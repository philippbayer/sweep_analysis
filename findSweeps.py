#!/usr/bin/env python3
"""
Selective Sweep Detection Pipeline for Selection Experiment (Optimized)
Identifies selective sweeps in Large and Small selection treatments vs Controls
"""

import subprocess
import os
import sys
from pathlib import Path
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

class SelectionExperimentAnalyzer:
    def __init__(self, vcf_dir=".", output_dir="sweep_analysis", window_size=50000, 
                 step_size=10000, threads=None, ploidy=1):
        """
        Initialize the selection experiment sweep detection pipeline
        
        Args:
            vcf_dir: Directory containing VCF files
            output_dir: Output directory for results
            window_size: Window size for sliding window analysis (bp)
            step_size: Step size for sliding windows (bp)
            threads: Number of parallel threads (default: CPU count - 1)
            ploidy: Ploidy level (1 for haploid, 2 for diploid)
        """
        self.vcf_dir = Path(vcf_dir)
        self.output_dir = Path(output_dir)
        self.window_size = window_size
        self.step_size = step_size
        self.threads = threads or max(1, multiprocessing.cpu_count() - 1)
        self.ploidy = ploidy
        self.output_dir.mkdir(exist_ok=True)
        
        print(f"Analyzing {'HAPLOID' if ploidy == 1 else 'DIPLOID'} organism")
        
        self.vcf_files = self.organize_vcf_files()
        
    def organize_vcf_files(self):
        """Organize VCF files by lineage and treatment"""
        files = {'large': [], 'small': [], 'control': []}
        lineages = defaultdict(dict)
        
        for vcf in self.vcf_dir.glob("*.vcf*"):
            name = vcf.stem.replace('.vcf', '').replace('.gz', '')
            
            if '_L' in name:
                treatment = 'large'
                lineage = name.replace('_L', '')
            elif '_S' in name:
                treatment = 'small'
                lineage = name.replace('_S', '')
            else:
                treatment = 'control'
                lineage = name
            
            files[treatment].append(str(vcf))
            lineages[lineage][treatment] = str(vcf)
        
        print("Found VCF files:")
        print(f"  Large selection (_L): {len(files['large'])} files")
        print(f"  Small selection (_S): {len(files['small'])} files")
        print(f"  Control (no suffix): {len(files['control'])} files")
        print(f"  Total: {len(files['large']) + len(files['small']) + len(files['control'])} files")
        
        print(f"\nLineage breakdown:")
        for lineage in sorted(lineages.keys()):
            treatments = []
            if 'large' in lineages[lineage]:
                treatments.append('L')
            if 'small' in lineages[lineage]:
                treatments.append('S')
            if 'control' in lineages[lineage]:
                treatments.append('C')
            print(f"  {lineage}: {', '.join(treatments)}")
        
        print(f"\nUsing {self.threads} threads for parallel processing")
        
        self.lineages = lineages
        return files
    
    def check_dependencies(self):
        """Check if required tools are installed"""
        required_tools = ['vcftools', 'bcftools', 'bgzip', 'tabix']
        missing = []
        
        for tool in required_tools:
            if subprocess.run(['which', tool], capture_output=True).returncode != 0:
                missing.append(tool)
        
        if missing:
            print(f"ERROR: Missing required tools: {', '.join(missing)}")
            print("\nInstallation instructions:")
            print("  conda install -c bioconda vcftools bcftools htslib")
            return False
        return True
    
    def get_sample_ids(self, vcf_file):
        """Extract sample IDs from VCF file"""
        cmd = ['bcftools', 'query', '-l', vcf_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        samples = result.stdout.strip().split('\n')
        return samples
    
    def create_sample_file(self, vcf_files, output_file):
        """Create a file with sample IDs from multiple VCF files"""
        all_samples = set()
        for vcf in vcf_files:
            samples = self.get_sample_ids(vcf)
            all_samples.update(samples)
        
        with open(output_file, 'w') as f:
            for sample in sorted(all_samples):
                f.write(f"{sample}\n")
        
        return output_file
    
    def prepare_vcf(self, vcf_file):
        """Compress and index VCF if needed (run once per file)"""
        if not vcf_file.endswith('.gz'):
            vcf_gz = f"{vcf_file}.gz"
            if not Path(vcf_gz).exists():
                subprocess.run(['bgzip', '-c', vcf_file], 
                             stdout=open(vcf_gz, 'w'), check=True)
            vcf_file = vcf_gz
        
        if not Path(f"{vcf_file}.tbi").exists():
            subprocess.run(['tabix', '-p', 'vcf', vcf_file], check=True)
        
        return vcf_file
    
    def merge_vcf_files(self, vcf_files, output_file):
        """Merge multiple VCF files efficiently"""
        if len(vcf_files) == 1:
            return self.prepare_vcf(vcf_files[0])
        
        print(f"Merging {len(vcf_files)} VCF files...")
        
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            prepared = list(executor.map(self.prepare_vcf, vcf_files))
        
        cmd = ['bcftools', 'merge', '-o', output_file, '--threads', str(self.threads)] + prepared
        subprocess.run(cmd, check=True)
        
        return output_file
    
    def calculate_tajimas_d(self, vcf_file, label):
        """Calculate Tajima's D in sliding windows
        Note: For haploids, Tajima's D is less informative as there's no within-individual variation"""
        output_prefix = self.output_dir / f"{label}_tajima"
        
        cmd = [
            'vcftools',
            '--vcf', vcf_file,
            '--TajimaD', str(self.window_size),
            '--out', str(output_prefix)
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        return f"{output_prefix}.Tajima.D"
    
    def calculate_pi(self, vcf_file, label):
        """Calculate nucleotide diversity (pi) in sliding windows"""
        output_prefix = self.output_dir / f"{label}_pi"
        
        cmd = [
            'vcftools',
            '--vcf', vcf_file,
            '--window-pi', str(self.window_size),
            '--window-pi-step', str(self.step_size),
            '--out', str(output_prefix)
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        return f"{output_prefix}.windowed.pi"
    
    def calculate_fst(self, vcf_file, pop1_file, pop2_file, label):
        """Calculate Fst between populations"""
        output_prefix = self.output_dir / f"{label}_fst"
        
        cmd = [
            'vcftools',
            '--vcf', vcf_file,
            '--weir-fst-pop', pop1_file,
            '--weir-fst-pop', pop2_file,
            '--fst-window-size', str(self.window_size),
            '--fst-window-step', str(self.step_size),
            '--out', str(output_prefix)
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        return f"{output_prefix}.windowed.weir.fst"
    
    def calculate_ld(self, vcf_file, label, max_distance=50000):
        """Calculate linkage disequilibrium (r^2) within windows
        High LD suggests reduced recombination and potential selective sweeps"""
        output_prefix = self.output_dir / f"{label}_ld"
        
        cmd = [
            'vcftools',
            '--vcf', vcf_file,
            '--hap-r2',  # Use haplotype r^2 (appropriate for haploids)
            '--ld-window-bp', str(max_distance),
            '--out', str(output_prefix)
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        return f"{output_prefix}.hap.ld"
    
    def summarize_ld_by_window(self, ld_file, label):
        """Summarize LD statistics into genomic windows"""
        import pandas as pd
        import numpy as np
        
        print(f"  Summarizing LD for {label}...")
        
        try:
            ld_data = pd.read_csv(ld_file, sep='\t')
        except:
            print(f"  Warning: Could not read LD file {ld_file}")
            return None
        
        if len(ld_data) == 0:
            print(f"  Warning: No LD data in {ld_file}")
            return None
        
        windows = []
        
        chroms = ld_data['CHR'].unique()
        
        for chrom in chroms:
            chrom_data = ld_data[ld_data['CHR'] == chrom]
            
            max_pos = max(chrom_data['POS1'].max(), chrom_data['POS2'].max())
            
            for window_start in range(0, int(max_pos), self.step_size):
                window_end = window_start + self.window_size
                
                window_ld = chrom_data[
                    ((chrom_data['POS1'] >= window_start) & (chrom_data['POS1'] < window_end)) |
                    ((chrom_data['POS2'] >= window_start) & (chrom_data['POS2'] < window_end))
                ]
                
                if len(window_ld) > 0:
                    windows.append({
                        'CHROM': chrom,
                        'BIN_START': window_start,
                        'BIN_END': window_end,
                        'POS': window_start + self.window_size // 2,
                        'MEAN_R2': window_ld['R^2'].mean(),
                        'MEDIAN_R2': window_ld['R^2'].median(),
                        'MAX_R2': window_ld['R^2'].max(),
                        'N_PAIRS': len(window_ld)
                    })
        
        summary = pd.DataFrame(windows)
        
        output_file = self.output_dir / f"{label}_ld_summary.txt"
        summary.to_csv(output_file, sep='\t', index=False)
        
        print(f"LD summary saved: {output_file}")
        return output_file
    
    def calculate_statistics_parallel(self, vcf_file, label):
        """Calculate Tajima's D, pi, and LD in parallel"""
        with ProcessPoolExecutor(max_workers=3) as executor:
            tajima_future = executor.submit(self.calculate_tajimas_d, vcf_file, label)
            pi_future = executor.submit(self.calculate_pi, vcf_file, label)
            ld_future = executor.submit(self.calculate_ld, vcf_file, label)
            
            tajima_file = tajima_future.result()
            pi_file = pi_future.result()
            ld_file = ld_future.result()
            
            ld_summary = self.summarize_ld_by_window(ld_file, label)
        
        return tajima_file, pi_file, ld_summary
    
    def analyze_comparison(self, treatment1_vcfs, treatment2_vcfs, comparison_name):
        """Analyze a specific treatment comparison"""
        print(f"\n{'='*60}")
        print(f"ANALYZING: {comparison_name}")
        print(f"{'='*60}")
        
        results_dir = self.output_dir / comparison_name.replace(' ', '_').replace('vs', 'vs')
        results_dir.mkdir(exist_ok=True)
        
        print("Step 1/4: Merging VCF files...")
        merged1 = results_dir / "treatment1_merged.vcf"
        merged2 = results_dir / "treatment2_merged.vcf"
        
        with ProcessPoolExecutor(max_workers=2) as executor:
            future1 = executor.submit(self.merge_vcf_files, treatment1_vcfs, str(merged1))
            future2 = executor.submit(self.merge_vcf_files, treatment2_vcfs, str(merged2))
            
            merged1 = future1.result()
            merged2 = future2.result()
        
        print("Step 2/4: Creating sample files...")
        pop1_samples = results_dir / "pop1_samples.txt"
        pop2_samples = results_dir / "pop2_samples.txt"
        
        self.create_sample_file([merged1], str(pop1_samples))
        self.create_sample_file([merged2], str(pop2_samples))
        
        print("Step 3/4: Calculating statistics...")
        combined_vcf = results_dir / "combined.vcf"
        combined_vcf = self.merge_vcf_files([merged1, merged2], str(combined_vcf))
        
        with ProcessPoolExecutor(max_workers=3) as executor:
            # Tajima, pi, and LD for treatment 1
            stats1_future = executor.submit(
                self.calculate_statistics_parallel, 
                merged1, 
                f"{comparison_name.split()[0]}"
            )
            
            # Tajima, pi, and LD for treatment 2
            stats2_future = executor.submit(
                self.calculate_statistics_parallel,
                merged2,
                f"{comparison_name.split()[2]}"
            )
            
            # Fst
            fst_future = executor.submit(
                self.calculate_fst,
                str(combined_vcf),
                str(pop1_samples),
                str(pop2_samples),
                comparison_name.replace(' ', '_')
            )
            
            tajima1, pi1, ld1 = stats1_future.result()
            tajima2, pi2, ld2 = stats2_future.result()
            fst = fst_future.result()
        
        print("Step 4/4: Identifying candidate regions...")
        self.identify_sweep_candidates(tajima1, tajima2, pi1, pi2, ld1, ld2, 
                                      fst, comparison_name, results_dir)
        
        return results_dir
    
    def identify_sweep_candidates(self, tajima1, tajima2, pi1, pi2, fst_file, 
                                  comparison_name, output_dir):
        """Identify candidate sweep regions for this comparison"""
        import pandas as pd
        import numpy as np
        
        taj1 = pd.read_csv(tajima1, sep='\t')
        taj2 = pd.read_csv(tajima2, sep='\t')
        pi1_data = pd.read_csv(pi1, sep='\t')
        pi2_data = pd.read_csv(pi2, sep='\t')
        fst = pd.read_csv(fst_file, sep='\t')
        
        taj1 = taj1[taj1['TajimaD'].notna()]
        taj2 = taj2[taj2['TajimaD'].notna()]
        pi1_data = pi1_data[pi1_data['PI'].notna()]
        pi2_data = pi2_data[pi2_data['PI'].notna()]
        fst = fst[fst['WEIGHTED_FST'].notna()]
        
        taj1['POS'] = (taj1['BIN_START'] + self.window_size // 2).astype(int)
        taj2['POS'] = (taj2['BIN_START'] + self.window_size // 2).astype(int)
        pi1_data['POS'] = ((pi1_data['BIN_START'] + pi1_data['BIN_END']) // 2).astype(int)
        pi2_data['POS'] = ((pi2_data['BIN_START'] + pi2_data['BIN_END']) // 2).astype(int)
        fst['POS'] = ((fst['BIN_START'] + fst['BIN_END']) // 2).astype(int)
        
        merged = (taj1[['CHROM', 'POS', 'TajimaD']].rename(columns={'TajimaD': 'TajimaD_1'})
                  .merge(taj2[['CHROM', 'POS', 'TajimaD']].rename(columns={'TajimaD': 'TajimaD_2'}), 
                         on=['CHROM', 'POS'], how='outer')
                  .merge(pi1_data[['CHROM', 'POS', 'PI']].rename(columns={'PI': 'PI_1'}), 
                         on=['CHROM', 'POS'], how='outer')
                  .merge(pi2_data[['CHROM', 'POS', 'PI']].rename(columns={'PI': 'PI_2'}), 
                         on=['CHROM', 'POS'], how='outer')
                  .merge(fst[['CHROM', 'POS', 'WEIGHTED_FST']], 
                         on=['CHROM', 'POS'], how='outer'))
        
        merged['Delta_TajimaD'] = merged['TajimaD_1'] - merged['TajimaD_2']
        merged['Delta_PI'] = merged['PI_1'] - merged['PI_2']
        
        candidates = merged.copy()
        
        fst_threshold = np.percentile(candidates['WEIGHTED_FST'].dropna(), 95) if candidates['WEIGHTED_FST'].notna().sum() > 0 else 0
        taj1_threshold = np.percentile(candidates['TajimaD_1'].dropna(), 5) if candidates['TajimaD_1'].notna().sum() > 0 else 0
        pi1_threshold = np.percentile(candidates['PI_1'].dropna(), 5) if candidates['PI_1'].notna().sum() > 0 else 0
        
        high_fst = candidates['WEIGHTED_FST'] > fst_threshold
        low_taj1 = candidates['TajimaD_1'] < taj1_threshold
        low_pi1 = candidates['PI_1'] < pi1_threshold
        
        strong_candidates = candidates[high_fst & (low_taj1 | low_pi1)].copy()
        strong_candidates = strong_candidates.sort_values('WEIGHTED_FST', ascending=False)
        
        output_file = output_dir / "candidate_sweep_regions.txt"
        with open(output_file, 'w') as f:
            f.write(f"# Candidate Selective Sweep Regions\n")
            f.write(f"# Comparison: {comparison_name}\n")
            f.write(f"# Window size: {self.window_size} bp\n")
            f.write(f"# Criteria: High Fst (top 5%) AND (Low Tajima's D OR Low pi in treatment 1)\n")
            f.write(f"#\n")
            f.write(f"# Found {len(strong_candidates)} candidate regions\n\n")
            
            if len(strong_candidates) > 0:
                f.write(strong_candidates.to_string(index=False))
            else:
                f.write("No strong candidates found with current thresholds.\n")
        
        merged_file = output_dir / "all_statistics.csv.gz"
        merged.to_csv(merged_file, index=False, compression='gzip')
        
        print(f"Found {len(strong_candidates)} strong candidate regions")
        
        return strong_candidates
    
    def run_all_comparisons(self):
        """Run all possible comparisons"""
        
        if not self.check_dependencies():
            return
        
        print("="*60)
        print("SELECTION EXPERIMENT SWEEP ANALYSIS")
        print("="*60)
        
        comparisons = []
        
        # Large vs Control
        if self.vcf_files['large'] and self.vcf_files['control']:
            comparisons.append(('Large vs Control', 
                              self.vcf_files['large'], 
                              self.vcf_files['control']))
        
        # Small vs Control
        if self.vcf_files['small'] and self.vcf_files['control']:
            comparisons.append(('Small vs Control', 
                              self.vcf_files['small'], 
                              self.vcf_files['control']))
        
        # Large vs Small
        if self.vcf_files['large'] and self.vcf_files['small']:
            comparisons.append(('Large vs Small', 
                              self.vcf_files['large'], 
                              self.vcf_files['small']))
        
        results = {}
        
        for comparison_name, treatment1, treatment2 in comparisons:
            try:
                result_dir = self.analyze_comparison(treatment1, treatment2, comparison_name)
                results[comparison_name] = result_dir
            except Exception as e:
                print(f"ERROR in {comparison_name}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        self.generate_summary_report(results)
        
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE!")
        print(f"Results saved in: {self.output_dir}")
        print("="*60)
        
        return results
    
    def generate_summary_report(self, results):
        """Generate a summary report across all comparisons"""
        summary_file = self.output_dir / "SUMMARY_REPORT.txt"
        
        with open(summary_file, 'w') as f:
            f.write("="*60 + "\n")
            f.write("SELECTIVE SWEEP ANALYSIS SUMMARY\n")
            f.write("="*60 + "\n\n")
            
            f.write("EXPERIMENTAL DESIGN:\n")
            f.write(f"  Lineages analyzed: {', '.join(sorted(self.lineages.keys()))}\n")
            f.write(f"  Large selection replicates: {len(self.vcf_files['large'])}\n")
            f.write(f"  Small selection replicates: {len(self.vcf_files['small'])}\n")
            f.write(f"  Control replicates: {len(self.vcf_files['control'])}\n")
            f.write(f"  Window size: {self.window_size} bp\n")
            f.write(f"  Step size: {self.step_size} bp\n")
            f.write(f"  Threads used: {self.threads}\n\n")
            
            f.write("COMPARISONS ANALYZED:\n")
            for comparison_name, result_dir in results.items():
                f.write(f"\n{comparison_name}:\n")
                candidate_file = result_dir / "candidate_sweep_regions.txt"
                if candidate_file.exists():
                    with open(candidate_file) as cf:
                        lines = cf.readlines()
                        for line in lines:
                            if line.startswith('# Found'):
                                f.write(f"  {line.strip()}\n")
                                break
                f.write(f"  Results directory: {result_dir}\n")
            
            f.write("\n" + "="*60 + "\n")
            f.write("INTERPRETATION:\n")
            f.write("="*60 + "\n")
            
            if self.ploidy == 1:
                f.write("""
Selective sweep signatures in HAPLOID organisms:

1. High Fst between treatments - MOST IMPORTANT for haploids
   - Indicates divergent selection between treatments
   - Values > 0.15-0.2 suggest strong differential selection

2. Low nucleotide diversity (pi) - IMPORTANT
   - Indicates reduced genetic variation in swept regions
   - Compare to genome-wide baseline

3. High linkage disequilibrium (LD r^2) - IMPORTANT for sweep detection
   - High LD indicates reduced recombination or recent selective sweep
   - Selective sweeps create extended LD due to hitchhiking
   - Hard sweeps: Very high LD across region
   - Soft sweeps: Moderate LD with multiple haplotypes
   - Reference: "Consequences of recombination for the evolution of 
     the mating type locus in Chlamydomonas reinhardtii" shows LD
     predicts selection efficacy

4. Tajima's D - LESS INFORMATIVE for haploids
   - In haploids, all variation is between individuals (no heterozygosity)
   - Can still show deviation from neutral expectations
   - Negative values may indicate recent population growth or sweeps

KEY DIFFERENCES FOR HAPLOIDS:
- No within-individual genetic variation (no heterozygosity)
- All alleles are directly exposed to selection
- Selective sweeps can fix faster than in diploids
- Focus primarily on: Fst + Low pi + High LD

STRONGEST CANDIDATES show concordant signals:
- High Fst (divergent selection)
+ Low pi (reduced diversity) 
+ High LD (extended linkage / recent sweep)

These are "extra strong" candidates where multiple independent
signals point to the same genomic region under selection.

For Large vs Control: Regions under selection for large body size
For Small vs Control: Regions under selection for small body size  
For Large vs Small: Regions differentially selected between treatments

INTERPRETING LD PATTERNS:
- Very high LD (r^2 > 0.8) across >50kb: Classic hard sweep signature
- Moderate LD (r^2 0.4-0.7): Possible soft sweep or older sweep
- Low LD despite high Fst: Selection with maintained recombination
  (may indicate balancing selection or recent gene conversion)

Next steps:
- Examine genes in high-Fst + high-LD candidate regions (strongest evidence)
- Look for overlapping candidates across comparisons
- Consider allele frequency trajectories in these regions
- Functional annotation of candidate loci
- Test whether LD patterns suggest hard vs soft sweeps
""")
            else:
                f.write("""
Selective sweep signatures:
1. High Fst between treatments - indicates divergent selection
2. Low Tajima's D - indicates recent selective sweep
3. Low nucleotide diversity (pi) - indicates reduced diversity

Strong candidates show multiple concordant signals indicating
selection specifically acting in one treatment vs another.

For Large vs Control: Regions under selection for large body size
For Small vs Control: Regions under selection for small body size  
For Large vs Small: Regions differentially selected between treatments

Next steps:
- Examine genes in candidate regions
- Look for overlapping candidates across comparisons
- Validate with additional population genetic tests
- Consider functional annotation of candidate loci
""")
        
        print(f"\nSummary report saved to: {summary_file}")


if __name__ == "__main__":
    import sys
    
    vcf_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else None
    
    analyzer = SelectionExperimentAnalyzer(
        vcf_dir=vcf_dir,
        output_dir="sweep_analysis_results",
        window_size=50000,
        step_size=10000,
        threads=threads,
        ploidy=1             # 1 for haploid, 2 for diploid
    )
    
    results = analyzer.run_all_comparisons()
    
    print("\n" + "="*60)
    print("RESULTS SUMMARY:")
    print("="*60)
    for comparison, result_dir in results.items():
        print(f"{comparison}: {result_dir}")
    print("="*60)
