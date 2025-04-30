mkdir -p figs-for-paper

function copyFile() {
    cp figs/$1.pdf figs-for-paper/$2_$1.pdf
}

copyFile divide_and_conquer_schematic fig1
copyFile data_plots fig2
copyFile time_plot fig3
copyFile accuracy_plots fig4
copyFile scran_basics figS1
copyFile fixnu-chen figS2a
copyFile fixnu-diff-chen figS2b
copyFile fixnu-ibarra-soria-SM figS3a
copyFile fixnu-diff-ibarra-soria-SM figS3b
copyFile fixnu-ibarra-soria-PSM figS4a
copyFile fixnu-diff-ibarra-soria-PSM figS4b
copyFile buettner_elbo figS5a
copyFile chen_elbo figS5b
copyFile tung_elbo figS5c
copyFile zeisel_elbo figS5d
copyFile geweke_mu_all figS6
copyFile geweke_delta_all figS7
copyFile geweke_epsilon_all figS8
copyFile ess_mu_all figS9
copyFile ess_delta_all figS10
copyFile ess_epsilon_all figS11
copyFile cell_splitting figS12
copyFile point_estimates_mu figS13
copyFile point_estimates_delta figS14
copyFile point_estimates_epsilon figS15
copyFile norm_plot figS16
copyFile hpd_width_mu figS17
copyFile hpd_width_delta figS18
copyFile hpd_width_epsilon figS19

cp tables/data-summary.tex figs-for-paper/tab1_data-summary.tex
cp tables/hmc-comparison.tex figs-for-paper/tab2_hmc-comparison.tex
