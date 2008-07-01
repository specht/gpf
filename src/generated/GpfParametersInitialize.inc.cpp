// GPF parameter intializations
mk_Parameters[r_GpfParameterName::Masses] = new k_GpfParameterEnum("masses", "Masses", "", r_MassesOption::Monoisotopic, gk_MassesOptions, gk_MassesOptionsReverse);
mk_Parameters[r_GpfParameterName::Protease] = new k_GpfParameterEnum("protease", "Protease", "", r_ProteaseOption::Trypsin, gk_ProteaseOptions, gk_ProteaseOptionsReverse);
mk_Parameters[r_GpfParameterName::MassError] = new k_GpfParameterInt("massError", "Precision (ppm)", "Specify the mass precision of the mass spectrometer here.", 700);
mk_Parameters[r_GpfParameterName::SearchSimilar] = new k_GpfParameterEnum("searchSimilar", "Similarity search", "Choose whether GPF should search for hits that are only partially equal to the query. Similar hits may be very different from the original query, but the mass is still correct.", r_YesNoOption::Yes, gk_YesNoOptions, gk_YesNoOptionsReverse);
mk_Parameters[r_GpfParameterName::SearchIntrons] = new k_GpfParameterEnum("searchIntrons", "Intron search", "", r_YesNoOption::Yes, gk_YesNoOptions, gk_YesNoOptionsReverse);
mk_Parameters[r_GpfParameterName::MaxIntronLength] = new k_GpfParameterInt("maxIntronLength", "Max intron length", "", 2100);
mk_Parameters[r_GpfParameterName::MinChainLength] = new k_GpfParameterInt("minChainLength", "Minimum aa chain length", "", 5);
mk_Parameters[r_GpfParameterName::FullDetails] = new k_GpfParameterEnum("fullDetails", "Full details", "Return full hit details", r_YesNoOption::No, gk_YesNoOptions, gk_YesNoOptionsReverse);

// fill GPF enum option arrays
gk_YesNoOptions["no"] = 0;
gk_YesNoOptionsReverse[0] = "no";
gk_YesNoOptions["yes"] = 1;
gk_YesNoOptionsReverse[1] = "yes";
gk_ProteaseOptions["trypsin"] = 0;
gk_ProteaseOptionsReverse[0] = "trypsin";
gk_MassesOptions["monoisotopic"] = 0;
gk_MassesOptionsReverse[0] = "monoisotopic";
gk_SearchIntronOptions["never"] = 0;
gk_SearchIntronOptionsReverse[0] = "never";
gk_SearchIntronOptions["always"] = 1;
gk_SearchIntronOptionsReverse[1] = "always";
gk_SearchIntronOptions["ifNecessary"] = 2;
gk_SearchIntronOptionsReverse[2] = "ifNecessary";
