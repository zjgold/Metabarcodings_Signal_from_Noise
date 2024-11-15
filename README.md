# Signal and noise in metabarcoding data

[Link to manuscript](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0285674)

Zachary Gold<sup>1,2</sup>, Andrew Olaf Shelton<sup>2</sup>, Helen R. Casendino<sup>3</sup>, Joe Duprey<sup>3</sup>, Ramón Gallego<sup>2</sup>, Amy Van Cise<sup>2</sup>, Mary Fisher<sup>4</sup>, Alexander J. Jensen<sup>2</sup>, Erin D’Agnese<sup>3</sup>, Elizabeth Andruszkiewicz Allan<sup>3</sup>, Ana Ramón-Laca<sup>2</sup>, Maya Garber-Yonts3, Michaela Labare<sup>5</sup>, Kim M. Parsons<sup>2</sup>, Ryan P. Kelly<sup>3</sup>

<sup>1</sup> Cooperative Institute for Climate, Ocean, & Ecosystem Studies, UW, Seattle, WA
<sup>2</sup> Northwest Fisheries Science Center, NMFS/NOAA, Seattle, WA
<sup>3</sup> School of Marine and Environmental Affairs, UW, Seattle, WA
<sup>4</sup> School of Aquatic Fisheries Science, UW, Seattle, WA
<sup>5</sup> Scripps Institution of Oceanography, UCSD, La Jolla


## Abstract
Metabarcoding is a powerful molecular tool for simultaneously surveying hundreds to thousands of species from a single sample, underpinning microbiome and environmental DNA (eDNA) methods. Deriving quantitative estimates of underlying biological communities from metabarcoding is critical for enhancing the utility of such approaches for health and conservation. Recent work has demonstrated that correcting for amplification biases in genetic metabarcoding data can yield quantitative estimates of template DNA concentrations. However, a major source of uncertainty in metabarcoding data stems from non-detections across technical PCR replicates where one replicate fails to detect a species observed in other replicates. Such non-detections are a special case of variability among technical replicates in metabarcoding data. While many sampling and amplification processes underlie observed variation in metabarcoding data, understanding the causes of non-detections is an important step in distinguishing signal from noise in metabarcoding studies. Here, we use both simulated and empirical data to 1) suggest how non-detections may arise in metabarcoding data, 2) outline steps to recognize uninformative data in practice, and 3) identify the conditions under which amplicon sequence data can reliably detect underlying biological signals. We show with both simulations and empirical data that, for a given species, the rate of non-detections among technical replicates is a function of both the template DNA concentration and species-specific amplification efficiency. Consequently, we conclude metabarcoding datasets are strongly affected by (1) deterministic amplification biases during PCR and (2) stochastic sampling of amplicons during sequencing—both of which we can model—but also by (3) stochastic sampling of rare molecules prior to PCR, which remains a frontier for quantitative metabarcoding. Our results highlight the importance of estimating species-specific amplification efficiencies and critically evaluating patterns of non-detection in metabarcoding datasets to better distinguish environmental signal from the noise inherent in molecular detections of rare targets.


## Description
This page is dedicated to hosting code generated for the Signal from Noise Manuscript currently in submission to PLOS Biology and will be made available as a pre-print.  Included on this page is
1. **Code**
    1. *calcofi_signal_noise_20220820.Rmd* This script does most of the analyses of empirical data sets and generates figures 2 and 3 in the paper.
    2. *mc31_organization_20210105.Rmd* This script organizes the mock community data.
    3. *mc31_organization_coastal_even_redo_20220408.Rmd* This script organizes additional mock community data.
    4. *taxonomy_matcher_12S_20210106.Rmd*  This script creates the final taxonomic paths for the mock community data.
    5. *taxonomy_matcher_12S_mock_even_redo_20220408.Rmd* This script creates the final taxonomic paths for the mock community data.
2. **Data**
    1. *All_amp_efficiencies-2022-06-03.csv* Calculated amplification efficiencies from the mock community data.
    2. *input_dna_conc_communities_20210103.csv* Mock community metadata including starting input concentrations of DNA.
    3. *microscopy_tech_nReads.RDS* Microscopy data from [Gold et al. 2022](https://github.com/zjgold/CalCOFI_eDNA). See manuscript for full description of data and how data were generated.
    4. *mifish_mock_community_data.RDS* Mock community data.
    5. *mifish_tech_nReads.RDS* Metabarcoding data from [Gold et al. 2022](https://github.com/zjgold/CalCOFI_eDNA). See manuscript for full description of data and how data were generated.

    *mock_sequences*
        1. *CRUX_DB*
            1. global and local reference databases from [Gold et al. 2021](https://github.com/zjgold/FishCARD)
            2. Output fasta files from taxonomy_matcher*.rmd scripts
            3. Blast output from salmon sequences.
        2. *hash.key_updated_c19.RDS* Final updated taxonomy table after resolving conflicts.
        3. *hash.key_updated.RDS* Anacapa derived taxonomy table.
        4. *mock_even_redo*
            1.*c19_fishcard_ASV_raw_taxonomy_60.txt* Anacapa output for mock community using the global reference database
            2. *metadata_kenai1_20220408.csv* metadata file
        5. *updated_0106*
            1. *12S_fishcard_taxonomy_tables* Anacapa output for mock community using the local Calfiornia Current Large Marine Ecosystem reference database
            2. *c19_fishcard_taxonomy_tables* Anacapa output for mock community using the Global reference database
            3. *p16S_shark_taxonomy_tables* Anacapa output for mock community using a 16S fish reference database (data not used in this manuscript)


Github will be updated with pre-print, NCBI SRA, and Dryad information as they are generated and made available through the review process.
