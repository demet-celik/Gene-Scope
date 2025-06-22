# Gene Scope

**Gene Scope** is a user-friendly bioinformatics tool for performing a complete DNA-to-protein analysis pipeline. Features include RNA secondary structure prediction, multiple sequence alignment with ClustalW, and parsed BLAST results using real biological data. Implemented as a Jupyter Notebook (`.ipynb`), this project guides users through the essential steps of DNA and protein analysis, from raw sequence to evolutionary insights. The notebook is designed to be run on **Google Colab** for maximum accessibility and ease of use.

---

## Project Overview

Gene Scope is structured to walk users through the essential steps of bioinformatics analysis, from raw DNA to protein function and evolutionary relationships. The project is divided into two main parts:

### 1. DNA Sequence Analysis

#### a. Reading FASTA Files
- **Function:** `read_fasta`
- **Concept:** FASTA is a standard text-based format for representing nucleotide or peptide sequences. Reading FASTA files is the first step in almost any bioinformatics workflow.
- **Project example files:**  
  For this project, example FASTA files have already been downloaded and are included in a separate folder within the GitHub project repository. These files will be used as sample data for the analyses in Gene Scope.  
  **Before starting the project, users must download these example FASTA files from the GitHub repository and upload them to their Google Colab session when running the `Gene_Scope.ipynb` notebook.**
- **Where to find FASTA files:**  
  The example FASTA files used in this project—such as human genome DNA sequences and IFT140 protein sequences from human, mouse, and zebrafish—are real biological data that researchers use in their work. These files can also be downloaded from major biological databases, including:
  - **NCBI (National Center for Biotechnology Information):**  
    [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)  
    Search for your gene or protein of interest, navigate to the sequence record, and use the “Send to” or “Download” options to obtain the FASTA file.
  - **Ensembl Genome Browser:**  
    [https://www.ensembl.org/](https://www.ensembl.org/)  
    Provides comprehensive genome data for a wide range of species. Use the search bar to find your gene/protein, then download the sequence in FASTA format.
  - **UniProt:**  
    [https://www.uniprot.org/](https://www.uniprot.org/)  
    A leading resource for protein sequence and functional information. Search for your protein and download the FASTA sequence from the entry page.
- **How researchers use these platforms:**  
  These platforms are essential resources for both biologists and bioinformaticians:
  - **Biologists** use these databases to retrieve gene and protein sequences for experimental design, comparative studies, and functional annotation.
  - **Bioinformaticians** rely on these resources for large-scale data mining, sequence analysis, evolutionary studies, and to develop computational tools and pipelines.
- **Why included:** To teach users how to handle real biological data and prepare it for analysis.

#### b. Nucleotide Composition & GC Content
- **Function:** `dna_analysis`
- **Concept:**  
  - **Nucleotide counts** (A, T, C, G) provide basic information about the sequence.
  - **GC content** is the percentage of guanine (G) and cytosine (C) bases in a DNA or RNA molecule. It is calculated as:  
    ```
    GC content (%) = ((Number of G + Number of C) / Total number of bases) × 100
    ```
- **Why GC content matters (real-world usages):**
  - **Thermal Stability:** DNA with higher GC content is more thermally stable because G-C pairs form three hydrogen bonds (compared to two for A-T pairs). This property is crucial for organisms living in high-temperature environments (thermophiles), whose genomes often have elevated GC content.
  - **PCR Primer Design:** In molecular biology labs, GC content is a key factor in designing primers for PCR (Polymerase Chain Reaction). Primers with very high or very low GC content can lead to inefficient or non-specific amplification.
  - **Genome Annotation:** GC content is used to identify gene-rich regions, as coding sequences in many organisms tend to have higher GC content than non-coding regions. For example, in the human genome, GC-rich isochores often correspond to gene-dense areas.
  - **Species Identification and Evolution:** Different species and even different genomic regions within a species can have characteristic GC content. This can be used for taxonomic classification, evolutionary studies, and to detect horizontal gene transfer events.
  - **Sequencing Quality Control:** In next-generation sequencing, abnormal GC content can indicate contamination or sequencing bias, so it is routinely checked in quality control pipelines.
- **Why included:**  
  GC content is a fundamental metric in genomics and molecular biology, with applications ranging from basic research to clinical diagnostics. Including GC content analysis in this project helps users understand its broad significance and practical applications in real-world biological research.

#### c. Transcription (DNA to mRNA)
- **Function:** `transcribe_dna`
- **Concept:**  
  Transcription is the process by which DNA is converted into messenger RNA (mRNA). In most bioinformatics workflows, the DNA sequence provided in a FASTA file is the **coding (sense) strand**. This makes transcription straightforward: you simply replace every thymine (T) with uracil (U) to get the mRNA sequence.

  **Why is this easy in most cases?**  
  - The coding (sense) strand has the same sequence as the mRNA (except T is replaced by U).
  - This is the convention for most public databases (NCBI, Ensembl, etc.), so you can usually transcribe directly.

  **Visual Chart: Coding vs. Template Strand**

  ```
  5' - Coding (Sense) Strand:      ATG GCT TGA ...
           | | | | | | | | | 
  3' - Template (Antisense) Strand: TAC CGA ACT ...
  ```

  - The **coding strand** (top) is what you usually get in a FASTA file.
  - The **template strand** (bottom) is used by RNA polymerase to synthesize mRNA, but is rarely provided directly.

  **Transcription Process:**
  - If you have the coding strand (most common):  
    - Replace T with U to get mRNA.
    - Example:  
      ```
      DNA (coding): 5' - ATG GCT TGA - 3'
      mRNA:         5' - AUG GCU UGA - 3'
      ```
  - If you have the template (antisense) strand (rare):  
    - First, find the **complementary sequence** (A↔T, C↔G) and reverse it to get the coding strand.
    - Then, transcribe as above.
    - Example:  
      ```
      DNA (template): 3' - TAC CGA ACT - 5'
      Complement:     5' - ATG GCT TGA - 3'  (coding strand)
      mRNA:           5' - AUG GCU UGA - 3'
      ```

  **Summary Table:**

  | Strand Type         | Sequence Example      | What To Do for mRNA         |
  |---------------------|----------------------|-----------------------------|
  | Coding (Sense)      | 5' - ATG GCT TGA - 3'| Replace T with U            |
  | Template (Antisense)| 3' - TAC CGA ACT - 5'| Find complement, then T→U   |

- **Why included:**  
  This step demonstrates the central dogma of molecular biology and prepares the sequence for further analysis. Understanding the difference between coding and template strands is crucial for accurate transcription and downstream analyses. In rare cases where a FASTA file contains the template strand, users must generate the complementary (coding) strand before transcribing.

#### d. mRNA Analysis (Codon Usage, Start/Stop Codons, Structure)
- **Function:** `mrna_analysis`
- **Concepts:**
  - **Codon usage:** The frequency of each codon in mRNA, which can affect translation efficiency and protein expression.
  - **Start/stop codons:** Indicate where translation begins and ends.
  - **Secondary structure prediction:**  
    - **What is mRNA secondary structure?**  
      mRNA is not just a linear string of nucleotides; it can fold back on itself to form local structures such as hairpins, loops, and stems due to base pairing (A-U and G-C). These structures are called **secondary structures**.
    - **Why does it matter?**  
      - Secondary structure can affect how efficiently an mRNA is translated into protein.
      - It can influence mRNA stability (how long it lasts in the cell).
      - Certain structures can regulate gene expression by hiding or exposing ribosome binding sites or regulatory elements.
      - Some viruses and regulatory RNAs rely on specific secondary structures for their function.
    - **How is it predicted?**  
      - Computational tools like **RNAfold** use thermodynamic models to predict the most likely folding pattern of an mRNA sequence.
      - The output is usually a “dot-bracket” notation and a minimum free energy (MFE) value, which estimates the stability of the structure.
    - **Example of dot-bracket notation:**  
      ```
      Sequence:      AUGGCUACGUGA
      Structure:     (((....)))..
      ```
      - Here, `(` and `)` indicate paired bases (forming a stem), and `.` indicates unpaired bases (loops or bulges).
    - **Simple ASCII diagram:**
      ```
      5' - A U G G C U A C G U G A - 3'
             | | |         | |
             C C G         G C
      ```
      - The vertical lines show base pairs forming a stem, while the unpaired bases form a loop.
    - **Interpreting the results:**  
      - The **minimum free energy (MFE)** value (in kcal/mol) tells you how stable the predicted structure is. More negative values mean a more stable structure.
      - The predicted structure can be visualized as a 2D diagram, helping researchers identify important features like hairpins or open regions.
    - **Real-world applications:**  
      - **mRNA vaccine design:**  
        In the development of mRNA vaccines (such as those for COVID-19), scientists carefully design the mRNA sequence not only to encode the correct protein (like the viral spike protein), but also to optimize its secondary structure.  
        - **Why?** The secondary structure affects how long the mRNA persists in the body, how efficiently it is translated into protein, and how the immune system recognizes it.
        - **How?** Researchers use computational tools to predict and modify the mRNA’s folding, avoiding strong hairpins or structures that could block ribosome access, and ensuring the mRNA is stable enough to survive in the cell but not so stable that it resists translation.
        - **Result:** Optimized mRNA structure leads to higher protein production and a stronger, more reliable immune response.
      - Studying regulatory elements in untranslated regions (UTRs) of mRNAs.
      - Understanding how mutations affect mRNA folding and function.

- **Why included:**  
  Predicting and understanding mRNA secondary structure is crucial for interpreting gene expression, designing experiments, and developing RNA-based therapeutics. Including this analysis in the project helps users appreciate the complexity and importance of RNA structure in molecular biology.

#### e. Translation (mRNA to Protein)
- **Function:** `translate_mrna`
- **Concept:**  
  Translation is the process by which the genetic code carried by mRNA is decoded to produce a specific sequence of amino acids, resulting in a functional protein. This process occurs in the ribosome, a molecular machine found in all living cells.

  **How translation works:**
  1. **Initiation:**  
     - The ribosome scans the mRNA for the first start codon (**AUG**), which codes for the amino acid methionine.
     - Translation always begins at the first AUG in the correct reading frame.
  2. **Elongation:**  
     - The ribosome reads the mRNA in sets of three nucleotides (codons).
     - Each codon specifies a particular amino acid, which is added to the growing polypeptide chain.
     - Transfer RNAs (tRNAs) bring the correct amino acids to the ribosome according to the codon sequence.
  3. **Termination:**  
     - When the ribosome encounters a stop codon (**UAA, UAG, or UGA**), translation ends.
     - The completed protein is released from the ribosome.

  **Reading frames and accuracy:**  
  - The mRNA can be read in three possible reading frames, but only one will produce the correct protein.
  - If translation starts at the wrong position, the resulting protein will be incorrect (frameshift).
  - The function in this project finds the first AUG and translates from there, stopping at the first in-frame stop codon.

  **Example:**
  ```
  mRNA: 5' - GGA AUG GCU UGA CCA - 3'
                |   |   |   |
              (start)  A   Stop
  Protein:     Met  Ala
  ```
  - Translation starts at the first AUG (Met), continues with Ala, and stops at UGA.

  **Why is this important?**
  - The correct translation of mRNA is essential for producing functional proteins.
  - Mutations that add or remove nucleotides can cause frameshifts, leading to nonfunctional or harmful proteins.
  - Many genetic diseases are caused by errors in translation or mutations that introduce premature stop codons.

  **Real-world applications:**
  - **Biotechnology:** Synthetic genes are designed with optimized start/stop codons for efficient protein production.
  - **Medicine:** Understanding translation is key to developing gene therapies and treating diseases caused by translation errors.
  - **Research:** Protein expression studies rely on accurate translation to study protein function and structure.

- **Why included:**  
  This step illustrates how genetic information is converted into functional proteins, a central concept in molecular biology. It also demonstrates the importance of reading frames, start/stop codons, and the consequences of mutations on protein synthesis.

#### f. Protein Analysis (Molecular Weight)
- **Function:** `protein_analysis`
- **Concept:**  
  The molecular weight (MW) of a protein is the sum of the atomic masses of all its amino acids, usually measured in Daltons (Da) or kilodaltons (kDa). Knowing the molecular weight is essential for many laboratory techniques and for understanding protein function.

- **Why is molecular weight important?**
  - **Laboratory techniques:**  
    - **SDS-PAGE (Sodium Dodecyl Sulfate Polyacrylamide Gel Electrophoresis):**  
      This is a common method used to separate proteins based on their size. By comparing the migration of a protein sample to a set of molecular weight standards, researchers can estimate the size of the protein they have purified or expressed.
    - **Western blotting:**  
      After SDS-PAGE, proteins are transferred to a membrane and detected with antibodies. Knowing the expected molecular weight helps confirm the identity of the protein band.
  - **Clinical diagnostics:**  
    - Some diseases are diagnosed by detecting abnormal proteins or fragments with altered molecular weights (e.g., truncated proteins in genetic disorders).
  - **Protein purification:**  
    - Molecular weight is used to select appropriate filters or columns for isolating proteins of interest.

- **Real-life example:**  
  - **Insulin production:**  
    Recombinant human insulin, used to treat diabetes, is produced in bacteria or yeast. After expression, the protein is purified and analyzed by SDS-PAGE. The expected molecular weight of insulin is about 5.8 kDa. If the band on the gel matches this size, it confirms successful production and purity of the insulin protein before it is formulated for medical use.

- **Why included:**  
  Calculating molecular weight is a fundamental step in protein analysis, bridging computational predictions with real-world laboratory and clinical applications.

---

### 2. Protein Sequence Analysis & Multiple Sequence Alignment (MSA)

#### a. Combining Multiple FASTA Files
- **Function:** `concatenate_fastas`
- **Concept:**  
  Combining sequences from different species (orthologs) allows for comparative analysis and evolutionary studies. In this project, you will work with three example FASTA files, each containing the IFT140 protein sequence from a different species: human, mouse, and zebrafish.

- **Note:**  
  **You must use the example FASTA files provided in the GitHub repository.** These files are named to indicate their species of origin (e.g., `IFT140_HUMAN_NCBI.fasta`, `IFT140_House_Mouse_NCBI.fasta`, `IFT140_Zebrafish_NCBI.fasta`).  
  Before running this part of the notebook, download these files from the GitHub repo and upload them to your Google Colab session.

- **About IFT140:**  
  IFT140 (Intraflagellar Transport 140) is a protein that plays a crucial role in the assembly and maintenance of cilia and flagella in eukaryotic cells. Mutations in IFT140 are associated with several human diseases, including ciliopathies.  
  These sequences were chosen for this project because:
  - IFT140 is highly conserved across species, making it ideal for demonstrating multiple sequence alignment (MSA) and evolutionary analysis.
  - It is biologically significant and relevant to both basic research and clinical studies.

- **Purpose in the workflow:**  
  - The three IFT140 FASTA files are concatenated into a single multi-FASTA file.
  - This combined file is then used for **multiple sequence alignment (MSA)** to identify conserved regions and evolutionary relationships.
  - After MSA, the user can select any of the three sequences (e.g., human, mouse, or zebrafish IFT140) for further analysis, such as running a BLAST search to find similar proteins in databases.

- **Why included:**  
  This step prepares real, biologically meaningful data for downstream analyses, demonstrating how comparative genomics and evolutionary studies are performed in bioinformatics.

#### b. Multiple Sequence Alignment (MSA)
- **Function:** `run_msa_with_clustalw`
- **Concept:** MSA aligns three or more protein sequences to identify regions of similarity, which may indicate functional, structural, or evolutionary relationships.
- **Why included:** MSA is fundamental for studying protein families, evolution, and functional domains.

#### c. Highlighting Conserved Regions and Interpreting MSA Results
- **Function:** `display_highlighted_msa`
- **Concept:**  
  After performing a multiple sequence alignment (MSA), it is important to identify and visualize **conserved regions**—positions in the alignment where the same amino acid (or nucleotide) is present across all sequences. These regions are often functionally or structurally important, as they have been preserved by evolution due to selective pressure. Interpreting the results of a multiple sequence alignment goes beyond just visualizing the alignment. This function quantifies and explains the biological significance of the observed conservation.

- **How it works in this project:**  
  - The function scans each column of the alignment and checks if all sequences have the same amino acid at that position (excluding gaps).
  - The indices of these conserved positions are collected.
  - The alignment is then displayed with a color-coded highlight for each conserved amino acid, making it easy to visually spot these regions.
  - A summary of conserved positions is also printed in a readable format.
  - Calculates the **number and percentage of conserved positions** in the alignment.
  - Provides a summary of the alignment, including the number of sequences, alignment length, and overall percent identity.
  - Offers a biological interpretation based on the degree of conservation:
    - **High conservation** (e.g., >70%): Indicates strong evolutionary pressure to maintain function/structure, suggesting the proteins are closely related and likely share critical roles.
    - **Moderate conservation** (e.g., 40–70%): Suggests shared ancestry and functional domains, but with some divergence.
    - **Low conservation** (<40%): Indicates more distant evolutionary relationships or functional divergence.
  - Explains the possible biological roles of conserved regions, such as:
    - Active sites of enzymes
    - Structural elements necessary for protein folding
    - Binding sites for other molecules (proteins, DNA, ligands)

- **Why is this important?**
  - **Functional significance:** Conserved regions often correspond to active sites, binding domains, or structural motifs essential for the protein’s function.
  - **Evolutionary insight:** High conservation suggests evolutionary constraints, indicating that changes in these regions may be deleterious.
  - **Research applications:** Identifying conserved regions helps in designing experiments (e.g., mutagenesis), drug targeting, and understanding disease mechanisms.
  - **Evolutionary biology:** Helps trace the evolutionary history of genes/proteins and understand how function is preserved or diversified.
  - **Biomedical research:** Identifies potential targets for drug design or genetic engineering.
  - **Functional annotation:** Guides experimental studies by highlighting regions likely to be important for activity or interaction.

- **Visualization example:**
  ```
  Sequence1:  M A V Y F D H R A E A P D S S G V P V L I S W H S S V C V L A V G S V N P S T G G C V D ...
  Sequence2:  M A L Y F D H R I K A P D T P S S P S H I T W H P T H P F L A V A S I S P S S G G N V D ...
  Sequence3:  M A V Y F D H R A E A P D S S G V P V L I S W H S S V C V L A V G S V N P S T G G C V D ...
                * *   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
  ```
  (Asterisks or color highlights indicate conserved positions.)

  - **Example output for the interpretation:**
  ```
  ▶ Alignment Overview:
    - Number of Sequences: 3
    - Total Alignment Length: 1481 residues

  ▶ Conservation Analysis:
    - Number of Conserved Positions: 320
    - Overall Percent Identity: 21.6%

  ▶ Biological Significance:
    Conserved regions are powerful indicators of biological importance. Because these
    amino acids have resisted change over evolutionary time, they often correspond to:
      - The active site of an enzyme.
      - Key structural components required for the protein to fold correctly.
      - Sites for binding to other molecules (like DNA or other proteins).

  ▶ Conclusion:
    This is a moderate degree of conservation, suggesting a shared ancestry and the presence of common functional or structural domains, though significant divergence has also occurred.
  ```

- **Why included:**  
  This step provides an intuitive, visual way to interpret the results of MSA, helping users quickly identify biologically meaningful patterns in the data and also translates raw alignment data into meaningful biological insights, making the results actionable and understandable for both beginners and experienced researchers.

#### d. BLAST Search
- **Function:** `blast_and_parse_sequence`
- **Concept:**  
  **BLAST (Basic Local Alignment Search Tool)** is one of the most widely used algorithms in bioinformatics for comparing an input sequence (query) against a database of known sequences. It identifies regions of local similarity, helping researchers find homologous genes or proteins, infer function, and explore evolutionary relationships.

- **What does BLAST do?**
  - **Sequence comparison:** BLAST rapidly aligns your query sequence to millions of sequences in public databases (such as NCBI’s non-redundant protein or nucleotide databases).
  - **Finds homologs:** It detects similar sequences (homologs) in other organisms, which can suggest shared ancestry or function.
  - **Functional annotation:** By finding matches to well-characterized proteins, BLAST helps predict the function of unknown sequences.
  - **Evolutionary analysis:** BLAST results can reveal how conserved a sequence is across species, and help build phylogenetic trees.
  - **Medical and research applications:** Used for pathogen identification, gene discovery, primer design, and more.

- **How it works in this project:**  
  - After MSA, the user can select one of the IFT140 protein sequences (human, mouse, or zebrafish) and use it as a query for BLAST.
  - The function runs a BLAST search, retrieves the top hits, and provides a detailed interpretation of the results, including percent identity, E-value, and functional annotation.

- **Why included:**  
  BLAST is a cornerstone of modern bioinformatics, enabling researchers to connect their findings to the vast body of biological knowledge in public databases.

---

## Why Google Colab?

- **Automatic setup:** Installs all required tools and packages.
- **No local configuration needed:** Just upload your FASTA files and run the notebook.
- **Free and accessible:** Anyone with a Google account can use it.

---

## Educational Value

Gene Scope is designed to:
- Introduce the logic and workflow of bioinformatics, from raw sequence data to functional and evolutionary insights.
- Explain the biological significance behind each computational step, including real-world applications such as mRNA vaccine design, protein analysis in the lab, and evolutionary studies.
- Provide hands-on experience with real biological data and widely used tools (e.g., MSA, BLAST, RNAfold).
- Help users develop practical skills for both research and applied bioinformatics, bridging the gap between computation and biology.

---

## Future Directions

- Integration of more advanced analyses (e.g., motif finding, phylogenetics).
- Support for additional file formats and more complex datasets.
- Enhanced visualizations and interactive exploration.
- Deeper biological interpretation and literature integration.
- **Potential to be developed into a web application** for broader accessibility and ease of use.

---
