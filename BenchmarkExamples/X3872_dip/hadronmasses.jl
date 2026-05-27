# input hadron masses
print("Mass units in MeV")
const unit_choice = 1000; 
const mπ0 = 0.134977unit_choice; const mπc = 0.13957061unit_choice; const mπ = (2.0*mπc + mπ0)/3; 
const meta = 0.547862unit_choice
const mkc = 0.493677unit_choice; const mk0 = 0.497611unit_choice; const mk = (mkc+mk0)/2
const mkstarc = 0.89176unit_choice; const mkstar0 = 0.89555unit_choice; const mkstar = (mkstarc+mkstar0)/2; 
Γkstar = 0.05unit_choice;

const mdc = 1.86965unit_choice; const md0 = 1.86483unit_choice; const md = (mdc+md0)/2
const mdstarc = 2.01026unit_choice; const mdstar0 = 2.00685unit_choice; const mdstar = (mdstarc+mdstar0)/2;

const md10 = 2.4208unit_choice; const md1c = 2.4232unit_choice; const md1 = (md10+md1c)/2; 
const md2 = 2.4607unit_choice; const md2c = 2.4654unit_choice; 
const Γd1=0.0317unit_choice; const Γd2 = 0.047unit_choice

const mds = 1.96834unit_choice; const mdsstar = 2.1122unit_choice
const mds1 = 2.5351unit_choice; const Γds1 = 0.00092unit_choice; 
const mds2 = 2.5691unit_choice; const Γds2 = 0.0169unit_choice

const mjψ = 3.0969unit_choice; const mηc = 2.9839unit_choice; 
const mχc0 = 3.41471unit_choice; const mχc1 = 3.51067unit_choice; const mhc = 3.52538unit_choice; mχc2 = 3.55617unit_choice
const mψp = 3.6861unit_choice; const mηcp = 3.6375unit_choice; 
const mψ_3770 = 3.7737unit_choice; const mψ2 = 3.8222unit_choice
const mψ3 = 3.84271unit_choice; 
const mχc1_3872 = 3.87169unit_choice; const mzc_3900 = 3.8884unit_choice
const mx_3915 = 3.9184unit_choice; const mχc2_3930 = 3.9222unit_choice; const mzc_4020 = 4.0241unit_choice
const mψ_4040 = 4.039unit_choice; const mχc1_4140 = 4.1468unit_choice; const mψ_4160 = 4.191unit_choice; 
const mψ_4230 = 4.22unit_choice
const mχc1_4274 = 4.274unit_choice; const mψ_4360 = 4.368unit_choice; mψ_4415 = 4.421unit_choice
const mzc_4430 = 4.478unit_choice; const mψ_4660 = 4.633unit_choice

const mηb = 9.3987unit_choice; const mΥ = 9.4603unit_choice; 
const mχb0 = 9.85944unit_choice
const mχb1 = 9.89278unit_choice; const mhb = 9.8993unit_choice; const mχb2 = 9.91221unit_choice; 
const mΥ_2s = 10.02326unit_choice
const mΥ2 = 10.1637unit_choice; 
const mχb0_2p = 10.2325unit_choice; const mχb1_2p = 10.25546unit_choice; const mχb2_2p = 10.26865unit_choice
const mΥ3 = 10.3552unit_choice; const mχb1_3p = 10.5134unit_choice; const mχb2_3p = 10.524unit_choice; 
const mΥ_4s = 10.5794unit_choice; const mzb_10610 = 10.6072unit_choice; const mzb_10650 = 10.6522unit_choice; 
const mΥ_10860 = 10.8852unit_choice; const mΥ_11020 = 11.0unit_choice


const mp = 0.938272081unit_choice; const mn = 0.939565413unit_choice;

const mΛ = 1.115683unit_choice; 
const mΣplus = 1.18937unit_choice; const mΣ0 = 1.192642unit_choice; const mΣminus = 1.197449unit_choice

const mΞ0 = 1.31486unit_choice; mΞminus = 1.32171unit_choice

const mΩ = 1.67245unit_choice

const mΛc = 2.28646unit_choice;
const mΣcplusplus = 2.45397unit_choice; const mΣcplus = 2.4529unit_choice; const mΣc0 = 2.45375unit_choice
const mΣcplusplus_2520 = 2.51841unit_choice; const mΣcplus_2520 = 2.5175unit_choice; const mΣc0_2520 = 2.51848unit_choice
const mΞcplus = 2.46794unit_choice; const mΞc0 = 2.47090unit_choice; 
const mΞcplusp = 2.5784unit_choice; const mΞc0p = 2.5792unit_choice
const mΩc = 2.6952unit_choice

const mΞccplusplus = 3.6212unit_choice

const mΛb = 5.6196unit_choice; 
const mΣbplus = 5.81056unit_choice; const mΣbminus = 5.81564unit_choice
const mΞb0 = 5.797unit_choice; const mΞbminus = 5.7919unit_choice
const mΩb = 6.0461unit_choice

const ħc = 0.197327unit_choice;

qcm(m,m1,m2) = sqrt((m^2 - (m1+m2)^2) * (m^2-(m1-m2)^2) + 0im)/(2m)

nothing
# ε for Feynman prescription (notebook ϵ-prescription smearing)
const ϵ = 1e-6
