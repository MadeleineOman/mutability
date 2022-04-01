annotation_summarise = {"transcribed":["antisense_RNA","CDS","exon","gene","lnc_RNA","mRNA","ncRNA","miRNA","tRNA","vault_RNA",
                                       "V_gene_segment","Y_RNA","RNase_MRP_RNA","RNase_P_RNA","rRNA","scRNA","snoRNA","snRNA",
                                       "telomerase_RNA","transcript","primary_transcript","cDNA_match","CDS","exon","gene",
                                       "transcript","start_codon","stop_codon","stop_codon_redefined_as_selenocysteine"], 
                         "regulatory": ["insulator","enhancer","enhancer_blocking_element","promoter","regulatory_region","silencer",
                                        "TATA_box","GC_rich_promoter_region","CAAT_signal","locus_control_region"],
                         "protein_binding":["replication_regulatory_region"," replication_start_site","origin_of_replication",
                                            "meiotic_recombination_region","mitotic_recombination_region","recombination_feature",
                                            "non_allelic_homologous_recombination_region","matrix_attachment_site","protein_binding_site"],
                         "ignored":["direct_repeat","dispersed_repeat","tandem_repeat","repeat_instability_region","repeat_region",
                                    "match","D_loop","response_element","sequence_alteration","sequence_comparison","sequence_feature",
                                    "sequence_secondary_structure"," transcriptional_cis_regulatory_region","Genomic",
                                    "imprinting_control_region","pseudogene","region","biological_region","CAGE_cluster",
                                    "DNaseI_hypersensitive_site","epigenetically_modified_region","nucleotide_cleavage_site",
                                    "nucleotide_motif","centromere","conserved_region","microsatellite","minisatellite",
                                    "mobile_genetic_element","chromosome_breakpoint"],
                         "UTR3":["three_prime_UTR"],
                         "UTR5":["five_prime_UTR"]}

annotation_conversion = {}
for new_label,old_label_list in annotation_summarise.items(): 
    for old_label in old_label_list: 
        annotation_conversion[old_label] = new_label

def annotation_priorityLabel(annotation_list): 
    """list -->str 
    take in the list of annotaitons (if multiple at a site)
    output the one string that is the priority
    requires numpy as np"""
    annotation_prioritize = {"UTR3":1, "UTR5":2,"transcribed":3,"regulatory":4,"protein_binding":5}
    
    annotation_ints = []
    for label in annotation_list: 
        annotation_ints.append(annotation_prioritize[label])
    priority_int = np.min(annotation_ints)
    #convert back to label 
    for label,int_eq in annotation_prioritize.items(): 
        if int_eq == priority_int: 
            priority_label = label
    return priority_label