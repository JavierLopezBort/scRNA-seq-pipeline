
#########################################################3
# SAMPLEX
#########################################################3

if samplex == "Taelman":
    bp = f'{bf}GEO/human/GSE181558_Taelman_embryo_gonads'

elif samplex == "Irie":
    bp = f'{bf}GEO/human/GSE223036_Irie_DMRT1/scRNAseq'

elif samplex == "Seita":
    bp = f"{bf}GEO/marmoset/PRJNA902353_marmoset_pgclc"

elif samplex == "Hermann":
    bp = f"{bf}GEO/human/GSE109037_Hermann"

elif samplex == "Persio":
    bp = f"{bf}GEO/human/GSE153947_Persio"

elif samplex == "Murat":
    bp = f"{bf}GEO/human/PRJEB48242_Murat"

elif samplex == "Ivynatal":
    bp = f"{bf}ivn_scRNA"

elif samplex == "Ivynatal_fb":
    bp = f"{bf}ivn_scRNA"

elif samplex == "Murat_chimp":
    bp = f"{bf}GEO/chimpanzee/PRJEB48241_Murat"

elif samplex == "Huang":
    bp = f"{bf}GEO/water_buffalo/GSE190477_testes_scRNAseq"

elif samplex == "Sosa":
    bp = f"{bf}GEO/mouse/GSE239488_rOvary_Sosa_scRNAseq"

elif samplex == "Smela":
    bp = f"{bf}GEO/human/GSE268385_Merrick_DMI_scRNAseq"

elif samplex == "Murat_bonobo":
    bp = f"{bf}GEO/bonobo/ERP132587_Murat_Bonobo"

elif samplex == "Murat_gorilla":
    bp = f"{bf}GEO/gorilla/ERP132581_testes_murat_scRNAseq"

elif samplex == "Murat_macaque":
    bp = f"{bf}GEO/macaque/PRJEB48240_Murat"

elif samplex == "Murat_gibbon":
    bp = f"{bf}GEO/gibbon/ERR713_scRNAseq_gibbon"

elif samplex == "Murat_marmoset":
    bp = f"{bf}GEO/marmoset/PRJEB48246_Murat"

elif samplex == "Murat_mouse":
    bp = f"{bf}GEO/mouse/ERP132577_murat_testes_scRNAseq"

elif samplex == "Murat_oppossum":
    bp = f"{bf}GEO/oppossum/ERP132579_Murat_tests_scRNAseq"

elif samplex == "Murat_platypus":
    bp = f"{bf}GEO/platypus/ERR713_testes_scRNAseq"

elif samplex == "culm_annot":
    bp = f"{bf}GEO/culm_annot"

elif samplex == "Murat_chicken":
    bp = f"{bf}GEO/chicken/PRJEB48234_Murat"

elif samplex == "Saitou":
    bp = f"{bf}GEO/human/GSE231812_Saitou"

elif samplex == "Reich":
    bp = f"{bf}GEO/human/GSE32689_oocytes"

elif samplex == "Ge":
    bp = f"{bf}GEO/mouse/GSE128553_embryo_female"

#########################################################3
# SAMPLE DICT
#########################################################3

sample_dict_raw = dict()

if samplex == "gonads":

    sample_dict_raw["Mesonephros"] = [
        "1_Mesonephros_2DPGC_0H_GEX_202201988_S1_L001",
        "2_Mesonephros_2DPGC_0H_GEX_202201988_S2_L001"
        ] # The two Mesonephros are different samples, they should be separated
    
    sample_dict_raw["sample120H"] = [
    "sample120H_GEX_S2_L001",
    "sample120H_GEX_S2_L002",
    "sample120H_GEX_S2_L004"
    ]

    sample_dict_raw["sample48H"] = [
    "sample48H_GEX_S1_L001",
    "sample48H_GEX_S1_L002",
    "sample48H_GEX_S1_L004"
    ]

    sample_dict_raw["fgonad_1tr"] = [
        "G1"
    ]

    sample_dict_raw["fgonad_2tr"] = [
        "G2"
    ]

    sample_dict_raw["mgonad_1tr"] = [
        "G3"
    ]

    sample_dict_raw["mgonad_2tr"] = [
        "G4"
    ]

    sample_dict_raw["GEXa"] = [
    "GEXa_S1_L001",
    "GEXa_S1_L002",
    "GEXa_S1_L004"
    ]

    sample_dict_raw["GEXb"] = [
    "GEXb_S2_L001",
    "GEXb_S2_L002",
    "GEXb_S2_L004"
    ]

elif samplex == "Seita":
    sample_dict_raw["iPSC"] = [
        "02_Yasu_10X_cjC10_iPSC_S2_L001"
    ]

    sample_dict_raw["D2_PGCLC"] = [
        "03_Yasu_10X_cjC10_D2_PGCLC_S3_L001"
    ]

    sample_dict_raw["D6_PGCLC"] = [
        "04_Yasu_10X_cjC10_D6_PGCLC_S4_L001"
    ]

elif samplex == "Irie":
    sample_dict_raw["4i_ESC"] = [
        "ASsml67_10X_1_S1_L001",
        "ASsml67_10X_1_S1_L002"
    ]

    sample_dict_raw["NANOS3+PGCLC"] = [
        "ASsml67_10X_2_S2_L001",
        "ASsml67_10X_2_S2_L002"
    ]

    sample_dict_raw["DMRT1+PGCLC"] = [
        "ASsml67_10X_3_S3_L001",
        "ASsml67_10X_3_S3_L002"
    ]

    sample_dict_raw["d8_DAZL+PGCLC"] = [
        "ASsml67_10X_4_S4_L001",
        "ASsml67_10X_4_S4_L002"
    ]

elif samplex == "Hermann":
    sample_dict_raw["Sta-put_scytes"] = [
        "hA17-11_S1_L001",
        "hA17-11_S1_L002",
        "hA17-11_S1_L003",
        "hA17-11_S1_L004"
    ]

    sample_dict_raw["sdyst_sgenic_3"] = [
        "hA17-3_S1_L001"
    ]

    sample_dict_raw["sdyst_sgenic_4"] = [
        "hA17-4_S1_L001"
    ]

    sample_dict_raw["sdyst_sgenic_5"] = [
        "hA17-5_S1_L004"
    ]

    sample_dict_raw["sgonia_1"] = [
        "hSPG17-1_S1_L008"
    ]

    sample_dict_raw["sgonia_2"] = [
        "hSPG17-2_S1_L006"
    ]

elif samplex == "Persio":
    sample_dict_raw["testis_31y"] = [
        "N1_S2_L004"
    ]

    sample_dict_raw["testis_33y"] = [
        "N2_S4_L004"
    ]

    sample_dict_raw["testis_55y"] = [
        "N3_S5_L004"
    ]

    sample_dict_raw["testis_39y"] = [
        "X"
    ]

    sample_dict_raw["testis_25y"] = [
        "X"
    ]

    sample_dict_raw["testis_36y"] = [
        "X"
    ]

elif samplex == "Murat":
    # I think the alignment and fastq files correspondies are wrongly assigned
    sample_dict_raw["testis_32y_s1"] = [
        "SN007_S1_L002", 
        "SN007_S1_L003", 
        "SN007_S1_L004", 
        "SN007_S2_L001", 
        "SN007_S2_L002", 
        "SN007_S2_L003", 
        "SN007_S2_L004", 
        "SN007_S3_L001", 
        "SN007_S3_L002", 
        "SN007_S3_L003", 
        "SN007_S3_L004", 
        "SN007_S4_L001", 
        "SN007_S4_L002", 
        "SN007_S4_L003", 
        "SN007_S4_L004" 
    ]

    sample_dict_raw["testis_32y_s2"] = [
        "SN011_S1_L002", 
        "SN011_S1_L003", 
        "SN011_S1_L004", 
        "SN011_S2_L001", 
        "SN011_S2_L002", 
        "SN011_S2_L003", 
        "SN011_S2_L004", 
        "SN011_S3_L001", 
        "SN011_S3_L002", 
        "SN011_S3_L003", 
        "SN011_S3_L004", 
        "SN011_S4_L001", 
        "SN011_S4_L002", 
        "SN011_S4_L003", 
        "SN011_S4_L004" 
    ]

    sample_dict_raw["testis_32y_s3"] = [
        "SN052_Seq1_S17_L001",
        "SN052_Seq1_S17_L002",
        "SN052_Seq1_S17_L003",
        "SN052_Seq1_S17_L004",
        "SN052_Seq1_S18_L001",
        "SN052_Seq1_S18_L002",
        "SN052_Seq1_S18_L003",
        "SN052_Seq1_S18_L004",
        "SN052_Seq1_S19_L001",
        "SN052_Seq1_S19_L002",
        "SN052_Seq1_S19_L003",
        "SN052_Seq1_S19_L004",
        "SN052_Seq1_S20_L001",
        "SN052_Seq1_S20_L002",
        "SN052_Seq1_S20_L003",
        "SN052_Seq1_S20_L004",
        "SN052_Seq2_S17_L001",
        "SN052_Seq2_S17_L002",
        "SN052_Seq2_S17_L003",
        "SN052_Seq2_S17_L004",
        "SN052_Seq2_S18_L001",
        "SN052_Seq2_S18_L002",
        "SN052_Seq2_S18_L003",
        "SN052_Seq2_S18_L004",
        "SN052_Seq2_S19_L001",
        "SN052_Seq2_S19_L002",
        "SN052_Seq2_S19_L003",
        "SN052_Seq2_S19_L004",
        "SN052_Seq2_S20_L001",
        "SN052_Seq2_S20_L002",
        "SN052_Seq2_S20_L003",
        "SN052_Seq2_S20_L004",

    ]

    sample_dict_raw["testis_28y_s1"] = [
        "SN111_S5_L001",
        "SN111_S5_L002",
        "SN111_S5_L003",
        "SN111_S5_L004",
        "SN111_S6_L001",
        "SN111_S6_L002",
        "SN111_S6_L003",
        "SN111_S6_L004",
        "SN111_S7_L001",
        "SN111_S7_L002",
        "SN111_S7_L003",
        "SN111_S7_L004",
        "SN111_S8_L001",
        "SN111_S8_L002",
        "SN111_S8_L003",
        "SN111_S8_L004",
    ]

    sample_dict_raw["testis_28y_s2"] = [
        "SN142_Seq1_S21_L001",
        "SN142_Seq1_S21_L002",
        "SN142_Seq1_S21_L003",
        "SN142_Seq1_S21_L004",
        "SN142_Seq1_S22_L001",
        "SN142_Seq1_S22_L002",
        "SN142_Seq1_S22_L003",
        "SN142_Seq1_S22_L004",
        "SN142_Seq1_S23_L001",
        "SN142_Seq1_S23_L002",
        "SN142_Seq1_S23_L003",
        "SN142_Seq1_S23_L004",
        "SN142_Seq1_S24_L001",
        "SN142_Seq1_S24_L002",
        "SN142_Seq1_S24_L003",
        "SN142_Seq1_S24_L004",
        "SN142_Seq2_S1_L001",
        "SN142_Seq2_S1_L002",
        "SN142_Seq2_S1_L003",
        "SN142_Seq2_S1_L004",
        "SN142_Seq2_S2_L001",
        "SN142_Seq2_S2_L002",
        "SN142_Seq2_S2_L003",
        "SN142_Seq2_S2_L004",
        "SN142_Seq2_S3_L001",
        "SN142_Seq2_S3_L002",
        "SN142_Seq2_S3_L003",
        "SN142_Seq2_S3_L004",
        "SN142_Seq2_S4_L001",
        "SN142_Seq2_S4_L002",
        "SN142_Seq2_S4_L003",
        "SN142_Seq2_S4_L004",
    ]

    
elif samplex == "Murat_chimp":
    #testis_45y_ch
    sample_dict_raw["SN074"] = [
        "SN074_S1_L001",
        "SN074_S1_L002", 
        "SN074_S1_L003", 
        "SN074_S1_L004", 
        "SN074_S2_L001",
        "SN074_S2_L002",
        "SN074_S2_L003",
        "SN074_S2_L004",
        "SN074_S3_L001",
        "SN074_S3_L002",
        "SN074_S3_L003",
        "SN074_S3_L004",
        "SN074_S4_L001",
        "SN074_S4_L002",
        "SN074_S4_L003",
        "SN074_S4_L004",
    ]
    #testis_14y_ch
    sample_dict_raw["SN112"] = [
        "SN112_S17_L001",
        "SN112_S17_L002",
        "SN112_S17_L003",
        "SN112_S17_L004",
        "SN112_S18_L001",
        "SN112_S18_L002",
        "SN112_S18_L003",
        "SN112_S18_L004",
        "SN112_S19_L001",
        "SN112_S19_L002",
        "SN112_S19_L003",
        "SN112_S19_L004",
        "SN112_S20_L001",
        "SN112_S20_L002",
        "SN112_S20_L003",
        "SN112_S20_L004",
    ]
    #testis_21y_ch
    sample_dict_raw["SN193"] = [
        "SN193_S1_L001",
        "SN193_S1_L002",
        "SN193_S1_L003",
        "SN193_S1_L004",
        "SN193_S2_L001",
        "SN193_S2_L002",
        "SN193_S2_L003",
        "SN193_S2_L004",
        "SN193_S3_L001",
        "SN193_S3_L002",
        "SN193_S3_L003",
        "SN193_S3_L004",
        "SN193_S4_L001",
        "SN193_S4_L002",
        "SN193_S4_L003",
        "SN193_S4_L004",
    ]

elif samplex == "Ivynatal":
    sample_dict_raw["poolGEX"] = [
        "poolGEX_S1_L001"
    ]

elif samplex == "Ivynatal_fb":
    sample_dict_raw["poolGEX_fb"] = [
        "poolfeaturebarcode_S1_L001"
    ]

elif samplex == "Huang":
    sample_dict_raw["testis_24mo_bf_s12"] = [
        "peri-Puberty_H5GHFDSXY_S12_L001",
    ]
    sample_dict_raw["testis_24mo_bf_s6"] = [
        "peri-Puberty_H5GHFDSXY_S6_L001",
    ]
    sample_dict_raw["testis_24mo_bf_s79"] = [
        "peri-Puberty_H5GHFDSXY_S79_L001",
    ]
    sample_dict_raw["testis_24mo_bf_s81"] = [
        "peri-Puberty_H5GHFDSXY_S81_L001",
    ]
    sample_dict_raw["testis_3mo_bf_s19"] = [
        "pre-Puberty_H5GHFDSXY_S19_L001",
    ]
    sample_dict_raw["testis_3mo_bf_s42"] = [
        "pre-Puberty_H5GHFDSXY_S42_L001",
    ]
    sample_dict_raw["testis_3mo_bf_s73"] = [
        "pre-Puberty_H5GHFDSXY_S73_L001",
    ]
    sample_dict_raw["testis_3mo_bf_s88"] = [
        "pre-Puberty_H5GHFDSXY_S88_L001",
    ]
elif samplex == "Sosa":
    sample_dict_raw["rep1"] = [
        "D6_aggregate_1_S1_L001"
    ]
    sample_dict_raw["rep2"] = [
        "D6_aggregate_2_S6_L001",
        "D6_aggregate_2_seq2_S6_L001",
    ]
    sample_dict_raw["rep3"] = [
        "D6_aggregate_3_S7_L001",
        "D6_aggregate_3_seq2_S5_L001",
    ]
    sample_dict_raw["rep4"] = [
        "D6_aggregate_4_S13_L002",
        "D6_aggregate_4_seq2_S13_L003",
        "D6_aggregate_4_seq3_S13_L004",
        "D6_aggregate_4_seq4_S1_L002",
    ]

elif samplex == "Smela":
    sample_dict_raw["day1"] = [
        "24047FL-02-01-01_S0_L001",
    ]
    sample_dict_raw["day2"] = [
        "24047FL-02-02-01_S0_L001",
    ]
    sample_dict_raw["day3"] = [
        "24047FL-02-03-01_S0_L001",
    ]
    sample_dict_raw["day4"] = [
        "24047FL-02-04-01_S0_L001",
    ]
    sample_dict_raw["day5"] = [
        "24047FL-02-05-01_S0_L001",
    ]
    sample_dict_raw["day6"] = [
        "24047FL-02-06-01_S0_L001",
    ]
    sample_dict_raw["day7"] = [
        "24047FL-02-07-01_S0_L001",
    ]
    sample_dict_raw["day8"] = [
        "24047FL-02-08-01_S0_L001",
    ]

elif samplex == "Murat_bonobo":
    sample_dict_raw["SN219"] = [
        "SN219_Seq1_S13_L001",
        "SN219_Seq1_S13_L002",
        "SN219_Seq1_S13_L003",
        "SN219_Seq1_S13_L004",
        "SN219_Seq1_S14_L001",
        "SN219_Seq1_S14_L002",
        "SN219_Seq1_S14_L003",
        "SN219_Seq1_S14_L004",
        "SN219_Seq1_S15_L001",
        "SN219_Seq1_S15_L002",
        "SN219_Seq1_S15_L003",
        "SN219_Seq1_S15_L004",
        "SN219_Seq1_S16_L001",
        "SN219_Seq1_S16_L002",
        "SN219_Seq1_S16_L003",
        "SN219_Seq1_S16_L004",
        "SN219_Seq2_S1_L001",
        "SN219_Seq2_S1_L002",
        "SN219_Seq2_S1_L003",
        "SN219_Seq2_S1_L004",
        "SN219_Seq2_S2_L001",
        "SN219_Seq2_S2_L002",
        "SN219_Seq2_S2_L003",
        "SN219_Seq2_S2_L004",
        "SN219_Seq2_S3_L001",
        "SN219_Seq2_S3_L002",
        "SN219_Seq2_S3_L003",
        "SN219_Seq2_S3_L004",
        "SN219_Seq2_S4_L001",
        "SN219_Seq2_S4_L002",
        "SN219_Seq2_S4_L003",
        "SN219_Seq2_S4_L004",
    ]
    sample_dict_raw["SN224"] = [
        "SN224_Seq1_S5_L001",
        "SN224_Seq1_S5_L002",
        "SN224_Seq1_S5_L003",
        "SN224_Seq1_S5_L004",
        "SN224_Seq2_S2_L001",
        "SN224_Seq2_S2_L002",
        "SN224_Seq2_S2_L003",
        "SN224_Seq2_S2_L004",
    ]

elif samplex == "Murat_gorilla":
    sample_dict_raw["SN180"] = [
        "SN180_S5_L001",
        "SN180_S5_L002",
        "SN180_S5_L003",
        "SN180_S5_L004",
        "SN180_S6_L001",
        "SN180_S6_L002",
        "SN180_S6_L003",
        "SN180_S6_L004",
        "SN180_S7_L001",
        "SN180_S7_L002",
        "SN180_S7_L003",
        "SN180_S7_L004",
        "SN180_S8_L001",
        "SN180_S8_L002",
        "SN180_S8_L003",
        "SN180_S8_L004",

    ]

    sample_dict_raw["SN223"] = [
        "SN223_Seq1_S4_L001",
        "SN223_Seq1_S4_L002",
        "SN223_Seq1_S4_L003",
        "SN223_Seq1_S4_L004",
        "SN223_Seq2_S1_L001",
        "SN223_Seq2_S1_L002",
        "SN223_Seq2_S1_L003",
        "SN223_Seq2_S1_L004",
    ]

elif samplex == "Murat_macaque":
    sample_dict_raw["SN116"] = [
        "SN116_Seq1_S9_L001",
        "SN116_Seq1_S9_L002",
        "SN116_Seq1_S9_L003",
        "SN116_Seq1_S9_L004",
        "SN116_Seq1_S10_L001",
        "SN116_Seq1_S10_L002",
        "SN116_Seq1_S10_L003",
        "SN116_Seq1_S10_L004",
        "SN116_Seq1_S11_L001",
        "SN116_Seq1_S11_L002",
        "SN116_Seq1_S11_L003",
        "SN116_Seq1_S11_L004",
        "SN116_Seq1_S12_L001",
        "SN116_Seq1_S12_L002",
        "SN116_Seq1_S12_L003",
        "SN116_Seq1_S12_L004",
        "SN116_Seq2_S5_L001",
        "SN116_Seq2_S5_L002",
        "SN116_Seq2_S5_L003",
        "SN116_Seq2_S5_L004",
        "SN116_Seq2_S6_L001",
        "SN116_Seq2_S6_L002",
        "SN116_Seq2_S6_L003",
        "SN116_Seq2_S6_L004",
        "SN116_Seq2_S7_L001",
        "SN116_Seq2_S7_L002",
        "SN116_Seq2_S7_L003",
        "SN116_Seq2_S7_L004",
        "SN116_Seq2_S8_L001",
        "SN116_Seq2_S8_L002",
        "SN116_Seq2_S8_L003",
        "SN116_Seq2_S8_L004"
    ]
    sample_dict_raw["SN143"] = [
        "SN143_Seq1_S25_L001",
        "SN143_Seq1_S25_L002",
        "SN143_Seq1_S25_L003",
        "SN143_Seq1_S25_L004",
        "SN143_Seq1_S26_L001",
        "SN143_Seq1_S26_L002",
        "SN143_Seq1_S26_L003",
        "SN143_Seq1_S26_L004",
        "SN143_Seq1_S27_L001",
        "SN143_Seq1_S27_L002",
        "SN143_Seq1_S27_L003",
        "SN143_Seq1_S27_L004",
        "SN143_Seq1_S28_L001",
        "SN143_Seq1_S28_L002",
        "SN143_Seq1_S28_L003",
        "SN143_Seq1_S28_L004",
        "SN143_Seq2_S5_L001",
        "SN143_Seq2_S5_L002",
        "SN143_Seq2_S5_L003",
        "SN143_Seq2_S5_L004",
        "SN143_Seq2_S6_L001",
        "SN143_Seq2_S6_L002",
        "SN143_Seq2_S6_L003",
        "SN143_Seq2_S6_L004",
        "SN143_Seq2_S7_L001",
        "SN143_Seq2_S7_L002",
        "SN143_Seq2_S7_L003",
        "SN143_Seq2_S7_L004",
        "SN143_Seq2_S8_L001",
        "SN143_Seq2_S8_L002",
        "SN143_Seq2_S8_L003",
        "SN143_Seq2_S8_L004"
    ]

elif samplex == "Murat_gibbon":
    sample_dict_raw["SN181"] = [
        "SN181_S9_L001",
        "SN181_S9_L002",
        "SN181_S9_L003",
        "SN181_S9_L004",
        "SN181_S10_L001",
        "SN181_S10_L002",
        "SN181_S10_L003",
        "SN181_S10_L004",
        "SN181_S11_L001",
        "SN181_S11_L002",
        "SN181_S11_L003",
        "SN181_S11_L004",
        "SN181_S12_L001",
        "SN181_S12_L002",
        "SN181_S12_L003",
        "SN181_S12_L004",

    ]
    sample_dict_raw["SN194"] = [
        "SN194_S5_L001",
        "SN194_S5_L002",
        "SN194_S5_L003",
        "SN194_S5_L004",
        "SN194_S6_L001",
        "SN194_S6_L002",
        "SN194_S6_L003",
        "SN194_S6_L004",
        "SN194_S7_L001",
        "SN194_S7_L002",
        "SN194_S7_L003",
        "SN194_S7_L004",
        "SN194_S8_L001",
        "SN194_S8_L002",
        "SN194_S8_L003",
        "SN194_S8_L004",
    ]

elif samplex == "Murat_marmoset":
    sample_dict_raw["SN117"] = [
        "SN117_Seq1_S13_L001",
        "SN117_Seq1_S13_L002",
        "SN117_Seq1_S13_L003",
        "SN117_Seq1_S13_L004",
        "SN117_Seq1_S14_L001",
        "SN117_Seq1_S14_L002",
        "SN117_Seq1_S14_L003",
        "SN117_Seq1_S14_L004",
        "SN117_Seq1_S15_L001",
        "SN117_Seq1_S15_L002",
        "SN117_Seq1_S15_L003",
        "SN117_Seq1_S15_L004",
        "SN117_Seq1_S16_L001",
        "SN117_Seq1_S16_L002",
        "SN117_Seq1_S16_L003",
        "SN117_Seq1_S16_L004",
        "SN117_Seq2_S9_L001",
        "SN117_Seq2_S9_L002",
        "SN117_Seq2_S9_L003",
        "SN117_Seq2_S9_L004",
        "SN117_Seq2_S10_L001",
        "SN117_Seq2_S10_L002",
        "SN117_Seq2_S10_L003",
        "SN117_Seq2_S10_L004",
        "SN117_Seq2_S11_L001",
        "SN117_Seq2_S11_L002",
        "SN117_Seq2_S11_L003",
        "SN117_Seq2_S11_L004",
        "SN117_Seq2_S12_L001",
        "SN117_Seq2_S12_L002",
        "SN117_Seq2_S12_L003",
        "SN117_Seq2_S12_L004",
    ]
    sample_dict_raw["SN130"] = [
        "SN130_Seq1_S9_L001",
        "SN130_Seq1_S9_L002",
        "SN130_Seq1_S9_L003",
        "SN130_Seq1_S9_L004",
        "SN130_Seq1_S10_L001",
        "SN130_Seq1_S10_L002",
        "SN130_Seq1_S10_L003",
        "SN130_Seq1_S10_L004",
        "SN130_Seq1_S11_L001",
        "SN130_Seq1_S11_L002",
        "SN130_Seq1_S11_L003",
        "SN130_Seq1_S11_L004",
        "SN130_Seq1_S12_L001",
        "SN130_Seq1_S12_L002",
        "SN130_Seq1_S12_L003",
        "SN130_Seq1_S12_L004",
        "SN130_Seq2_S1_L001",
        "SN130_Seq2_S1_L002",
        "SN130_Seq2_S1_L003",
        "SN130_Seq2_S1_L004",
        "SN130_Seq2_S2_L001",
        "SN130_Seq2_S2_L002",
        "SN130_Seq2_S2_L003",
        "SN130_Seq2_S2_L004",
        "SN130_Seq2_S3_L001",
        "SN130_Seq2_S3_L002",
        "SN130_Seq2_S3_L003",
        "SN130_Seq2_S3_L004",
        "SN130_Seq2_S4_L001",
        "SN130_Seq2_S4_L002",
        "SN130_Seq2_S4_L003",
        "SN130_Seq2_S4_L004",
    ]

elif samplex == "Murat_mouse":
    sample_dict_raw["SN090"] = [
        "SN090_Seq1_S9_L001",
        "SN090_Seq1_S9_L002",
        "SN090_Seq1_S9_L003",
        "SN090_Seq1_S9_L004",
        "SN090_Seq1_S10_L001",
        "SN090_Seq1_S10_L002",
        "SN090_Seq1_S10_L003",
        "SN090_Seq1_S10_L004",
        "SN090_Seq1_S11_L001",
        "SN090_Seq1_S11_L002",
        "SN090_Seq1_S11_L003",
        "SN090_Seq1_S11_L004",
        "SN090_Seq1_S12_L001",
        "SN090_Seq1_S12_L002",
        "SN090_Seq1_S12_L003",
        "SN090_Seq1_S12_L004",
        "SN090_Seq2_S5_L001",
        "SN090_Seq2_S5_L002",
        "SN090_Seq2_S5_L003",
        "SN090_Seq2_S5_L004",
        "SN090_Seq2_S6_L001",
        "SN090_Seq2_S6_L002",
        "SN090_Seq2_S6_L003",
        "SN090_Seq2_S6_L004",
        "SN090_Seq2_S7_L001",
        "SN090_Seq2_S7_L002",
        "SN090_Seq2_S7_L003",
        "SN090_Seq2_S7_L004",
        "SN090_Seq2_S8_L001",
        "SN090_Seq2_S8_L002",
        "SN090_Seq2_S8_L003",
        "SN090_Seq2_S8_L004",
    ]
    sample_dict_raw["SN115"] = [
        "SN115_Seq1_S5_L001",
        "SN115_Seq1_S5_L002",
        "SN115_Seq1_S5_L003",
        "SN115_Seq1_S5_L004",
        "SN115_Seq1_S6_L001",
        "SN115_Seq1_S6_L002",
        "SN115_Seq1_S6_L003",
        "SN115_Seq1_S6_L004",
        "SN115_Seq1_S7_L001",
        "SN115_Seq1_S7_L002",
        "SN115_Seq1_S7_L003",
        "SN115_Seq1_S7_L004",
        "SN115_Seq1_S8_L001",
        "SN115_Seq1_S8_L002",
        "SN115_Seq1_S8_L003",
        "SN115_Seq1_S8_L004",
        "SN115_Seq2_S1_L001",
        "SN115_Seq2_S1_L002",
        "SN115_Seq2_S1_L003",
        "SN115_Seq2_S1_L004",
        "SN115_Seq2_S2_L001",
        "SN115_Seq2_S2_L002",
        "SN115_Seq2_S2_L003",
        "SN115_Seq2_S2_L004",
        "SN115_Seq2_S3_L001",
        "SN115_Seq2_S3_L002",
        "SN115_Seq2_S3_L003",
        "SN115_Seq2_S3_L004",
        "SN115_Seq2_S4_L001",
        "SN115_Seq2_S4_L002",
        "SN115_Seq2_S4_L003",
        "SN115_Seq2_S4_L004",
        ]

elif samplex == "Murat_oppossum":
    sample_dict_raw["SN067"] = [
        "SN067_S9_L001",
        "SN067_S9_L002",
        "SN067_S9_L003",
        "SN067_S9_L004",
        "SN067_S10_L001",
        "SN067_S10_L002",
        "SN067_S10_L003",
        "SN067_S10_L004",
        "SN067_S11_L001",
        "SN067_S11_L002",
        "SN067_S11_L003",
        "SN067_S11_L004",
        "SN067_S12_L001",
        "SN067_S12_L002",
        "SN067_S12_L003",
        "SN067_S12_L004",
    ]
    sample_dict_raw["SN071"] = [
        "SN071_S1_L001",
        "SN071_S1_L002",
        "SN071_S1_L003",
        "SN071_S1_L004",
        "SN071_S2_L001",
        "SN071_S2_L002",
        "SN071_S2_L003",
        "SN071_S2_L004",
        "SN071_S3_L001",
        "SN071_S3_L002",
        "SN071_S3_L003",
        "SN071_S3_L004",
        "SN071_S4_L001",
        "SN071_S4_L002",
        "SN071_S4_L003",
        "SN071_S4_L004",
    ]
    sample_dict_raw["SN277"] = [
        "SN277_Seq1_S5_L001",
        "SN277_Seq1_S5_L002",
        "SN277_Seq1_S5_L003",
        "SN277_Seq1_S5_L004",
        "SN277_Seq2_S5_L001",
        "SN277_Seq2_S5_L002",
        "SN277_Seq2_S5_L003",
        "SN277_Seq2_S5_L004",
        "SN277_Seq2_S6_L001",
        "SN277_Seq2_S6_L002",
        "SN277_Seq2_S6_L003",
        "SN277_Seq2_S6_L004",
        "SN277_Seq2_S7_L001",
        "SN277_Seq2_S7_L002",
        "SN277_Seq2_S7_L003",
        "SN277_Seq2_S7_L004",
        "SN277_Seq2_S8_L001",
        "SN277_Seq2_S8_L002",
        "SN277_Seq2_S8_L003",
        "SN277_Seq2_S8_L004",
    ]

elif samplex == "Murat_platypus":
    sample_dict_raw["SN253"] = [
        "SN253_Seq1_S5_L001",
        "SN253_Seq1_S5_L002",
        "SN253_Seq1_S5_L003",
        "SN253_Seq1_S5_L004",
        "SN253_Seq2_S5_L001",
        "SN253_Seq2_S5_L002",
        "SN253_Seq2_S5_L003",
        "SN253_Seq2_S5_L004",
        "SN253_Seq2_S6_L001",
        "SN253_Seq2_S6_L002",
        "SN253_Seq2_S6_L003",
        "SN253_Seq2_S6_L004",
        "SN253_Seq2_S7_L001",
        "SN253_Seq2_S7_L002",
        "SN253_Seq2_S7_L003",
        "SN253_Seq2_S7_L004",
        "SN253_Seq2_S8_L001",
        "SN253_Seq2_S8_L002",
        "SN253_Seq2_S8_L003",
        "SN253_Seq2_S8_L004",
    ]
    sample_dict_raw["SN260"] = [
        "SN260_Seq1_S4_L001",
        "SN260_Seq1_S4_L002",
        "SN260_Seq1_S4_L003",
        "SN260_Seq1_S4_L004",
        "SN260_Seq2_S1_L001",
        "SN260_Seq2_S1_L002",
        "SN260_Seq2_S1_L003",
        "SN260_Seq2_S1_L004",
    ]

elif samplex == "Murat_chicken":
    sample_dict_raw["SN264"] = [
        "SN264_Seq1_S3_L001",
        "SN264_Seq1_S3_L002",
        "SN264_Seq1_S3_L003",
        "SN264_Seq1_S3_L004",
        "SN264_Seq2_S1_L001",
        "SN264_Seq2_S1_L002",
        "SN264_Seq2_S1_L003",
        "SN264_Seq2_S1_L004",
        "SN264_Seq2_S2_L001",
        "SN264_Seq2_S2_L002",
        "SN264_Seq2_S2_L003",
        "SN264_Seq2_S2_L004",
        "SN264_Seq2_S3_L001",
        "SN264_Seq2_S3_L002",
        "SN264_Seq2_S3_L003",
        "SN264_Seq2_S3_L004",
        "SN264_Seq2_S4_L001",
        "SN264_Seq2_S4_L002",
        "SN264_Seq2_S4_L003",
        "SN264_Seq2_S4_L004",
    ]

    sample_dict_raw["SN265"] = [
        "SN265_Seq1_S4_L001",
        "SN265_Seq1_S4_L002",
        "SN265_Seq1_S4_L003",
        "SN265_Seq1_S4_L004",
        "SN265_Seq2_S5_L001",
        "SN265_Seq2_S5_L002",
        "SN265_Seq2_S5_L003",
        "SN265_Seq2_S5_L004",
        "SN265_Seq2_S6_L001",
        "SN265_Seq2_S6_L002",
        "SN265_Seq2_S6_L003",
        "SN265_Seq2_S6_L004",
        "SN265_Seq2_S7_L001",
        "SN265_Seq2_S7_L002",
        "SN265_Seq2_S7_L003",
        "SN265_Seq2_S7_L004",
        "SN265_Seq2_S8_L001",
        "SN265_Seq2_S8_L002",
        "SN265_Seq2_S8_L003",
        "SN265_Seq2_S8_L004",
    ]

elif samplex == "Saitou":
    sample_dict_raw["585B1_BTAG_WT"] = [
        "SRR24451203"
    ]

    sample_dict_raw["585B1_BTAG_TET1_KO1"] = [
        "SRR24451202"
    ]

    sample_dict_raw["585B1_BTAG_TET1_KO2"] = [
        "SRR24451201"
    ]

    sample_dict_raw["NCLCN_AGVT_c11"] = [
        "SRR24451200"
    ]

    sample_dict_raw["NCLCN_AGVT_c56"] = [
        "SRR24451199"
    ]

    sample_dict_raw["NCLCN_AGVT_c86"] = [
        "SRR24451198"
    ]

    sample_dict_raw["NCLCN_AGVT_c117"] = [
        "SRR24451197"
    ]

elif samplex == "Ge":
    sample_dict_raw["E11.5"] = [
        "GR11_S1_L001"
    ]
    sample_dict_raw["E12.5"] = [
        "GR12_S1_L001"
    ]
    sample_dict_raw["E13.5"] = [
        "GR13_S1_L001"
    ]
    sample_dict_raw["E14.5"] = [
        "GR14_S1_L001"
    ]

sample_trim_wildcard = [value for key, values in sample_dict_raw.items() for value in values]

if config["paired_end"]:
    sample_dict = {key: [f"{bp}/1_trimming{bm_tr}/{value}{read_preffix}{{read}}{read_suffix}{fastq_suffix}" for value in values]
    for key, values in sample_dict_raw.items()}
else:
    sample_dict = {key: [f"{bp}/1_trimming{bm_tr}/{value}{fastq_suffix}" for value in values]
    for key, values in sample_dict_raw.items()}


sns = list(sample_dict_raw.keys())

#########################################################3
# GENDER DICT
#########################################################3

if samplex == "Taelman":
    gender_dict = {
        "fgonad_1tr": "female",
        "fgonad_2tr": "female",
        "mgonad_1tr": "male",
        "mgonad_2tr": "male",
        "GEXa": "female",
        "GEXb": "female",
        "sample48H": "male&female",
        "sample120H": "female",
        "Mesonephros": "male&female"
    }

elif samplex == "Seita":
    gender_dict = {
        "iPSCs": "None",
        "PGCLCs_d2": "None",
        "PGCLCs_d6": "None"
    }

elif samplex == "Hermann":
    gender_dict = {
        "Sta-put_scytes": "male",
        "sdyst_sgenic_3": "male",
        "sdyst_sgenic_4": "male",
        "sdyst_sgenic_5": "male",
        "sgonia_1": "male",
        "sgonia_2": "male",
    }

elif samplex == "Irie":
    gender_dict = {
        "4i_ESCs": "None",
        "PGCLCs_N3+": "None",
        "PGCLCs_DM+": "None",
        "PGCLCs_DZ+": "None"
    }

elif samplex == "Persio":
    gender_dict = {
        "testis_31y": "male",
        "testis_33y": "male",
        "testis_55y": "male",
        "testis_39y": "male",
        "testis_25y": "male",
        "testis_36y": "male"
    }

elif samplex == "Murat":
    gender_dict = {
        "testis_28y_s1": "male",
        "testis_28y_s2": "male",
        "testis_32y_s1": "male",
        "testis_32y_s2": "male",
        "testis_32y_s3": "male"
    }

elif samplex == "Murat_chimp":
    gender_dict = {
        "SN074": "male",
        "SN112": "male",
        "SN193": "male"
    }

elif samplex == "Huang":
    gender_dict = {
        "testis_24mo_bf_s12": "male",
        "testis_24mo_bf_s6": "male",
        "testis_24mo_bf_s79": "male",
        "testis_24mo_bf_s81": "male",
        "testis_3mo_bf_s19": "male",
        "testis_3mo_bf_s42": "male",
        "testis_3mo_bf_s73": "male",
        "testis_3mo_bf_s88": "male",
    }
elif samplex == "Murat_bonobo":
    gender_dict = {
        "SN219": "male",
        "SN224": "male",
    }

elif samplex == "Sosa":
    gender_dict = {
        "rep1": "None",
        "rep2": "None",
        "rep3": "None",
        "rep4": "None",
    }

elif samplex == "Smela":
    gender_dict = {
        "day1": "None",
        "day2": "None",
        "day3": "None",
        "day4": "None",
        "day5": "None",
        "day6": "None",
        "day7": "None",
        "day8": "None"
    }
elif samplex == "Murat_gorilla":
    gender_dict = {
        "SN180": "male",
        "SN223": "male",
    }

elif samplex == "Murat_macaque":
    gender_dict = {
        "SN116": "male",
        "SN143": "male",
    }

elif samplex == "Murat_gibbon":
    gender_dict = {
        "SN181": "male",
        "SN194": "male",
    }

elif samplex == "Murat_marmoset":
    gender_dict = {
        "SN117": "male",
        "SN130": "male",
    }

elif samplex == "Murat_mouse":
    gender_dict = {
        "SN090": "male",
        "SN115": "male",
    }

elif samplex == "Murat_oppossum":
    gender_dict = {
        "SN067": "male",
        "SN071": "male",
        "SN277": "male",
    }

elif samplex == "Murat_platypus":
    gender_dict = {
        "SN253": "male",
        "SN260": "male",
    }

elif samplex == "Murat_chicken":
    gender_dict = {
        "SN264": "male",
        "SN265": "male",
    }

elif samplex == "Saitou":
    gender_dict = {
        "585B1_BTAG_WT": "None",
        "585B1_BTAG_TET1_KO1": "None",
        "585B1_BTAG_TET1_KO2": "None",
        "NCLCN_AGVT_c11": "None",
        "NCLCN_AGVT_c56": "None",
        "NCLCN_AGVT_c86": "None",
        "NCLCN_AGVT_c117": "None",
    }
elif samplex == "Ivynatal":
    gender_dict = {
        "HTO_C0251": "None",
        "HTO_C0252": "None",
        "HTO_C0253": "None",
        "HTO_C0254": "None",
        "HTO_C0255": "None",
        "HTO_C0256": "None"
    }

elif samplex == "Ge":
    gender_dict = {
        "E11.5": "female",
        "E12.5": "female",
        "E13.5": "female",
        "E14.5": "female",
    }

#########################################################3
# STUDY
#########################################################3

study = samplex

if samplex == "Murat_chimp" or samplex == "Murat_marmoset" or samplex == "Murat_macaque" or samplex == "Murat_chicken" or samplex == "Murat_mouse" or samplex == "Murat_bonobo" or samplex == "Murat_gorilla" or samplex == "Murat_gibbon" or samplex == "Murat_oppossum" or samplex == "Murat_platypus":
    study = "Murat"

#########################################################3
# CELL TYPE
#########################################################3

if samplex == "Taelman":
    cell_type_dict = {
        "fgonad_1tr": "embryo",
        "fgonad_2tr": "embryo",
        "mgonad_1tr": "embryo",
        "mgonad_2tr": "embryo",
        "GEXa": "embryo",
        "GEXb": "embryo",
        "sample48H": "embryo",
        "sample120H": "embryo",
        "Mesonephros": "embryo"
    }
elif samplex == "Seita":
    cell_type_dict = {
        "iPSC": "induced",
        "D2_PGCLC": "induced",
        "D6_PGCLC": "induced"
    }

elif samplex == "Hermann":
    cell_type_dict = {
        "Sta-put_scytes": "adult",
        "sdyst_sgenic_3": "adult",
        "sdyst_sgenic_4": "adult",
        "sdyst_sgenic_5": "adult",
        "sgonia_1": "adult",
        "sgonia_2": "adult",
    }

elif samplex == "Irie":
    cell_type_dict = {
        "4i_ESC": "induced",
        "NANOS3+PGCLC": "induced",
        "DMRT1+PGCLC": "induced",
        "d8_DAZL+PGCLC": "induced"
    }

elif samplex == "Persio":
    cell_type_dict = {
        "testis_31y": "adult",
        "testis_33y": "adult",
        "testis_55y": "adult",
        "testis_39y": "adult",
        "testis_25y": "adult",
        "testis_36y": "adult"
    }

elif samplex == "Murat":
    cell_type_dict = {
        "testis_28y_s1": "adult",
        "testis_28y_s2": "adult",
        "testis_32y_s1": "adult",
        "testis_32y_s2": "adult",
        "testis_32y_s3": "adult"
    }

elif samplex == "Ivynatal":
    cell_type_dict = {
        "HTO_C0251": "induced",
        "HTO_C0252": "induced",
        "HTO_C0253": "induced",
        "HTO_C0254": "induced",
        "HTO_C0255": "induced",
        "HTO_C0256": "induced"
    }

elif samplex == "Murat_chimp":
    cell_type_dict = {
        "SN074": "adult",
        "SN112": "adult",
        "SN193": "adult"
    }

elif samplex == "Huang":
    cell_type_dict = {
        "testis_24mo_bf_s12": "adult",
        "testis_24mo_bf_s6": "adult",
        "testis_24mo_bf_s79": "adult",
        "testis_24mo_bf_s81": "adult",
        "testis_3mo_bf_s19": "pre-pub",
        "testis_3mo_bf_s42": "pre-pub",
        "testis_3mo_bf_s73": "pre-pub",
        "testis_3mo_bf_s88": "pre-pub",
    }

elif samplex == "Murat_bonobo":
    cell_type_dict = {
        "SN219": "adult",
        "SN224": "adult",
    }

elif samplex == "Sosa":
    cell_type_dict = {
        "rep1": "induced",
        "rep2": "induced",
        "rep3": "induced",
        "rep4": "induced",
    }

elif samplex == "Smela":
    cell_type_dict = {
        "day1": "induced",
        "day2": "induced",
        "day3": "induced",
        "day4": "induced",
        "day5": "induced",
        "day6": "induced",
        "day7": "induced",
        "day8": "induced"
    }

elif samplex == "Murat_gorilla":
    cell_type_dict = {
        "SN180": "adult",
        "SN223": "adult",
    }

elif samplex == "Murat_macaque":
    cell_type_dict = {
        "SN116": "adult",
        "SN143": "adult",
    }

elif samplex == "Murat_gibbon":
    cell_type_dict = {
        "SN181": "adult",
        "SN194": "adult",
    }

elif samplex == "Murat_marmoset":
    cell_type_dict = {
        "SN117": "adult",
        "SN130": "adult",
    }

elif samplex == "Murat_mouse":
    cell_type_dict = {
        "SN090": "adult",
        "SN115": "adult"
    }

elif samplex == "Murat_oppossum":
    cell_type_dict = {
        "SN067": "adult",
        "SN071": "adult",
        "SN277": "adult",
    }

elif samplex == "Murat_platypus":
    cell_type_dict = {
        "SN253": "adult",
        "SN260": "adult",
    }

elif samplex == "Murat_chicken":
    cell_type_dict = {
        "SN264": "adult",
        "SN265": "adult",
    }

elif samplex == "Saitou":
    cell_type_dict = {
        "585B1_BTAG_WT": "induced",
        "585B1_BTAG_TET1_KO1": "induced",
        "585B1_BTAG_TET1_KO2": "induced",
        "NCLCN_AGVT_c11": "induced",
        "NCLCN_AGVT_c56": "induced",
        "NCLCN_AGVT_c86": "induced",
        "NCLCN_AGVT_c117": "induced",
    }

elif samplex == "Ge":
    cell_type_dict = {
        "E11.5": "embryo",
        "E12.5": "embryo",
        "E13.5": "embryo",
        "E14.5": "embryo",
    }

#########################################################3
# SPECIE
#########################################################3

if samplex == "Taelman":
    specie = "human"
elif samplex == "Seita":
    specie = "marmoset"
elif samplex == "Hermann":
    specie = "human"
elif samplex == "Irie":
    specie = "human"
elif samplex == "Persio":
    specie = "human"
elif samplex == "Murat":
    specie = "human"
elif samplex == "Murat_chimp":
    specie = "chimp"
elif samplex == "Huang":
    specie = "buffalo"
elif samplex == "Sosa":
    specie = "mouse"
elif samplex == "Smela":
    specie = "human"
elif samplex == "Murat_bonobo":
    specie = "bonobo"
elif samplex == "Murat_gorilla":
    specie = "gorilla"
elif samplex == "Murat_macaque":
    specie = "macaque"
elif samplex == "Murat_gibbon":
    specie = "gibbon"
elif samplex == "Murat_marmoset":
    specie = "marmoset"
elif samplex == "Murat_mouse":
    specie = "mouse"
elif samplex == "Murat_oppossum":
    specie = "oppossum"
elif samplex == "Murat_platypus":
    specie = "platypus"
elif samplex == "Murat_chicken":
    specie = "chicken"
elif samplex == "Saitou":
    specie = "human"
elif samplex == "Ivynatal":
    specie = "human"
elif samplex == "Ge":
    specie = "mouse"
    
