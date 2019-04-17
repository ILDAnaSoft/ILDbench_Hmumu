#!/bin/bash

source /afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/init_ilcsoft.sh
cd /afs/desy.de/user/s/skawada/sonas_work/stage_nnh250/analysis/output

rm *.root *~ LEFT/*.root RIGHT/*.root

sleep 5

. hadd.sh

sleep 5

root -l -b -q SkimVar.C+\(\"alldata1.root\",\"alldata1_select.root\",\"dataTree\",\"evt:processid:xsec:higgsdecay1:higgsdecay2:MC_mumu_costh:h_nt_m:h_nq:h_ne:h_nm:h_nt:h_nn:h_nboson_m:type:tau_m:forff1:forff2:ff1_costh:ff2_costh:ff_nt_m:forffff1:forffff2:forffff3:forffff4:ffff1_costh:ffff2_costh:ffff3_costh:ffff4_costh:ffff1_t_costh:ffff2_t_costh:ffff3_t_costh:ffff4_t_costh:massff:mass1:mass2:MC_mumupair_mass:MC_muminus_parent:MC_muplus_parent:nq:ne:nm:nt:nn:no:nt_m:npfos:ntrks:esum:costh_missmom:maxPt:Pt:Pt_all:mass_vis:mass_inv:principalthrust:costh_thrustaxis:n_muplus:n_muminus:muplus_E:muminus_E:muplus_costh:muminus_costh:muplus_Pt:muminus_Pt:muplus_charge_costh:muminus_charge_costh:muplus_azi:muminus_azi:muplus_d0:muplus_z0:muplus_phi:muplus_omega:muplus_tanlambda:muplus_d0sig:muplus_z0sig:muplus_Chi2:muplus_Ndf:muplus_Chi2Ndf:muplus_dEdx:muplus_dEdxError:muplus_RadiusOfInnermostHit:muminus_d0:muminus_z0:muminus_phi:muminus_omega:muminus_tanlambda:muminus_d0sig:muminus_z0sig:muminus_Chi2:muminus_Ndf:muminus_Chi2Ndf:muminus_dEdx:muminus_dEdxError:muminus_RadiusOfInnermostHit:sum_charge_costh:leadingmu_E:subleadingmu_E:leadingmu_costh:subleadingmu_costh:mumu_E:mumu_mass:mumu_Pt:mumu_costh:mumu_acop:sigma_mumu_mass:recoilmass:mumu_costh_tobeam:mumu_mom_mag:n_ISR:pseudomass:n_highPt:mc_muminus_PDG:mc_muplus_PDG:parent_mc_muminus_PDG:parent_mc_muplus_PDG:esum_mc:costh_missmom_mc:Pt_all_mc:overlay_E_mc:mumu_mass_mc:mumu_costh_mc:n_highPt_mc:mumu_mass_mc_smear_1_3:mumu_mass_mc_smear_5_4:mumu_mass_mc_smear_3_4:mumu_mass_mc_smear_2_4:mumu_mass_mc_smear_1_4:mumu_mass_mc_smear_5_5:mumu_mass_mc_smear_3_5:mumu_mass_mc_smear_2_5:mumu_mass_mc_smear_1_5:mumu_mass_mc_smear_5_6:mumu_mass_mc_smear_3_6:mumu_mass_mc_smear_2_6:mumu_mass_mc_smear_1_6\"\)

sleep 5

root -l -b -q SkimVar.C+\(\"alldata2.root\",\"alldata2_select.root\",\"dataTree\",\"evt:processid:xsec:higgsdecay1:higgsdecay2:MC_mumu_costh:h_nt_m:h_nq:h_ne:h_nm:h_nt:h_nn:h_nboson_m:type:tau_m:forff1:forff2:ff1_costh:ff2_costh:ff_nt_m:forffff1:forffff2:forffff3:forffff4:ffff1_costh:ffff2_costh:ffff3_costh:ffff4_costh:ffff1_t_costh:ffff2_t_costh:ffff3_t_costh:ffff4_t_costh:massff:mass1:mass2:MC_mumupair_mass:MC_muminus_parent:MC_muplus_parent:nq:ne:nm:nt:nn:no:nt_m:npfos:ntrks:esum:costh_missmom:maxPt:Pt:Pt_all:mass_vis:mass_inv:principalthrust:costh_thrustaxis:n_muplus:n_muminus:muplus_E:muminus_E:muplus_costh:muminus_costh:muplus_Pt:muminus_Pt:muplus_charge_costh:muminus_charge_costh:muplus_azi:muminus_azi:muplus_d0:muplus_z0:muplus_phi:muplus_omega:muplus_tanlambda:muplus_d0sig:muplus_z0sig:muplus_Chi2:muplus_Ndf:muplus_Chi2Ndf:muplus_dEdx:muplus_dEdxError:muplus_RadiusOfInnermostHit:muminus_d0:muminus_z0:muminus_phi:muminus_omega:muminus_tanlambda:muminus_d0sig:muminus_z0sig:muminus_Chi2:muminus_Ndf:muminus_Chi2Ndf:muminus_dEdx:muminus_dEdxError:muminus_RadiusOfInnermostHit:sum_charge_costh:leadingmu_E:subleadingmu_E:leadingmu_costh:subleadingmu_costh:mumu_E:mumu_mass:mumu_Pt:mumu_costh:mumu_acop:sigma_mumu_mass:recoilmass:mumu_costh_tobeam:mumu_mom_mag:n_ISR:pseudomass:n_highPt:mc_muminus_PDG:mc_muplus_PDG:parent_mc_muminus_PDG:parent_mc_muplus_PDG:esum_mc:costh_missmom_mc:Pt_all_mc:overlay_E_mc:mumu_mass_mc:mumu_costh_mc:n_highPt_mc:mumu_mass_mc_smear_1_3:mumu_mass_mc_smear_5_4:mumu_mass_mc_smear_3_4:mumu_mass_mc_smear_2_4:mumu_mass_mc_smear_1_4:mumu_mass_mc_smear_5_5:mumu_mass_mc_smear_3_5:mumu_mass_mc_smear_2_5:mumu_mass_mc_smear_1_5:mumu_mass_mc_smear_5_6:mumu_mass_mc_smear_3_6:mumu_mass_mc_smear_2_6:mumu_mass_mc_smear_1_6\"\)

sleep 5

. hadd2.sh

sleep 5

root -l -b <<EOF
.L proc.C+
process()
.q
EOF

sleep 5

cp alldata_select.root LEFT/.
cp alldata_select.root RIGHT/.
