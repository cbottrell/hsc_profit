CREATE TABLE IF NOT EXISTS Morphologies_ProFit_HSC_Sersic(

    DatabaseID VARCHAR(255),
    SnapNum INT(32),
    SubfindID INT(32),
    Camera VARCHAR(255),
    Band VARCHAR(255),
    ProcessFlag INT(32),
    
    Redshift FLOAT,
    ApparentMagnitude FLOAT,
    RightAscension FLOAT,
    Declination FLOAT,
    Sersic_xcen FLOAT,
    Sersic_ycen FLOAT,
    Sersic_mag FLOAT,
    Sersic_re FLOAT,
    Sersic_nser FLOAT,
    Sersic_axrat FLOAT,
    Sersic_ang FLOAT,
    Sersic_xcen_med FLOAT,
    Sersic_xcen_std FLOAT,
    Sersic_xcen_p84 FLOAT,
    Sersic_xcen_m16 FLOAT,
    Sersic_ycen_med FLOAT,
    Sersic_ycen_std FLOAT,
    Sersic_ycen_p84 FLOAT,
    Sersic_ycen_m16 FLOAT,
    Sersic_mag_med FLOAT,
    Sersic_mag_std FLOAT,
    Sersic_mag_p84 FLOAT,
    Sersic_mag_m16 FLOAT,
    Sersic_re_med FLOAT,
    Sersic_re_std FLOAT,
    Sersic_re_p84 FLOAT,
    Sersic_re_m16 FLOAT,
    Sersic_nser_med FLOAT,
    Sersic_nser_std FLOAT,
    Sersic_nser_p84 FLOAT,
    Sersic_nser_m16 FLOAT,
    Sersic_ang_med FLOAT,
    Sersic_ang_std FLOAT,
    Sersic_ang_p84 FLOAT,
    Sersic_ang_m16 FLOAT,
    Sersic_axrat_med FLOAT,
    Sersic_axrat_std FLOAT,
    Sersic_axrat_p84 FLOAT,
    Sersic_axrat_m16 FLOAT,
    ReducedChiSquared FLOAT,
    ProFound_segID INT(32),
    ProFound_uniqueID INT(32),
    ProFound_xcen FLOAT,
    ProFound_ycen FLOAT,
    ProFound_xmax FLOAT,
    ProFound_ymax FLOAT,
    ProFound_sep FLOAT,
    ProFound_flux FLOAT,
    ProFound_mag FLOAT,
    ProFound_flux_app FLOAT,
    ProFound_mag_app FLOAT,
    ProFound_cenfrac FLOAT,
    ProFound_N50 FLOAT,
    ProFound_N90 FLOAT,
    ProFound_N100 FLOAT,
    ProFound_R50 FLOAT,
    ProFound_R90 FLOAT,
    ProFound_R100 FLOAT,
    ProFound_SB_N50 FLOAT,
    ProFound_SB_N90 FLOAT,
    ProFound_SB_N100 FLOAT,
    ProFound_xsd FLOAT,
    ProFound_ysd FLOAT,
    ProFound_covxy FLOAT,
    ProFound_corxy FLOAT,
    ProFound_con FLOAT,
    ProFound_asymm FLOAT,
    ProFound_flux_reflect FLOAT,
    ProFound_mag_reflect FLOAT,
    ProFound_semimaj FLOAT,
    ProFound_semimin FLOAT,
    ProFound_axrat FLOAT,
    ProFound_ang FLOAT,
    ProFound_signif FLOAT,
    ProFound_FPlim FLOAT,
    ProFound_flux_err FLOAT,
    ProFound_mag_err FLOAT,
    ProFound_flux_err_sky FLOAT,
    ProFound_flux_err_skyRMS FLOAT,
    ProFound_flux_err_shot FLOAT,
    ProFound_flux_err_cor FLOAT,
    ProFound_cor_seg FLOAT,
    ProFound_sky_mean FLOAT,
    ProFound_sky_sum FLOAT,
    ProFound_skyRMS_mean FLOAT,
    ProFound_Nedge INT(32),
    ProFound_Nsky INT(32),
    ProFound_Nobject INT(32),
    ProFound_Nborder INT(32),
    ProFound_Nmask INT(32),
    ProFound_edge_frac FLOAT,
    ProFound_edge_excess FLOAT,
    ProFound_flag_border INT(32),
    ProFound_iter INT(32),
    ProFound_origfrac FLOAT,
    ProFound_Norig INT(32),
    ProFound_skyseg_mean FLOAT,
    ProFound_flag_keep INT(32),
    RPetro_Elliptical FLOAT,
    RPetro_Circular FLOAT,
    RMax_Elliptical FLOAT,
    RMax_Circular FLOAT,
    R80_Elliptical FLOAT,
    R80_Circular FLOAT,
    R50_Elliptical FLOAT,
    R50_Circular FLOAT,
    R20_Elliptical FLOAT,
    R20_Circular FLOAT,
    Asymmetry_xcen FLOAT,
    Asymmetry_ycen FLOAT,
    Asymmetry FLOAT,
    AsymmetryBkgDens FLOAT,
    RMSAsymmetrySquared FLOAT,
    OuterAsymmetry FLOAT,
    ShapeAsymmetry FLOAT,
    AsymmetryNoAperture_xcen FLOAT,
    AsymmetryNoAperture_ycen FLOAT,
    AsymmetryNoAperture FLOAT,
    Concentration_Elliptical FLOAT,
    Concentration_Circular FLOAT,
    Smoothness FLOAT,
    SmoothnessBkgDens FLOAT,
    Gini FLOAT,
    M20 FLOAT,
    SB1kpc FLOAT,
    ResidualAsymmetry_xcen FLOAT,
    ResidualAsymmetry_ycen FLOAT,
    ResidualAsymmetry FLOAT,
    ResidualAsymmetryNoAperture_xcen FLOAT,
    ResidualAsymmetryNoAperture_ycen FLOAT,
    ResidualAsymmetryNoAperture FLOAT,
    
    PRIMARY KEY (DatabaseID)
);
