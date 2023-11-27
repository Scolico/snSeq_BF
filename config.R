###################
###Set Work Path###
###################


workpath='/mnt/snRNA_seq'
datapath='/mnt/data/10x'

setwd(workpath)

###init_config###----
setClass('config', slots = list(result.path = 'character',
                                data.path = 'character',
                                code.path = 'character',
                                project = 'character',
                                pc.num = 'numeric',
                                seed.use = 'numeric',
                                resolution = 'vector',
                                nfeatures = 'numeric',
                                ref.marker = 'vector',
                                mt.genes = 'vector',
                                rb.genes = 'vector',
                                hb.genes = 'vector',
                                sex.genes = 'vector',
                                cc.genes = 'list'),
         
        prototype = list(
           result.path = './snrna/result',
           data.path = './snrna/data',
           code.path = './snrna/code',
           project = 'forebrain',
           pc.num = 50,
           seed.use = 2021,
           resolution = seq(0.1, 2, 0.1),
           nfeatures = 2000,
           ref.marker = c('C1qc', 'Cldn5', 'Slc1a3', 'Sox10', 'Mobp','Pdgfra','Slc6a13','Mecom',
                          'Vcan', 'Rbfox3', 'Gad1', 'Slc17a6', 'Sox2', 'Hopx', 'Dlx6os1', 'Tmem212'),
           rb.genes = c('Rpl7', 'Rpl31', 'Rpl37a', 'Rps6kc1', 'Rpl7a', 'Rpl12', 'Rpl35',
                        'Rps21', 'Rpl22l1', 'Rps3a1', 'Rps27', 'Rpl34', 'Rps20', 'Rps6', 'Rps8',
                        'Rps6ka1', 'Rpl11', 'Rpl22', 'Rpl9', 'Rpl5', 'Rplp0', 'Rpl6', 'Rpl21',
                        'Rpl32', 'Rps9', 'Rpl28', 'Rps5', 'Rps19', 'Rps16', 'Rps11', 'Rpl13a',
                        'Rpl18', 'Rps17', 'Rps3', 'Rpl27a', 'Rps13', 'Rps15a', 'Rplp2',
                        'Rpl18a', 'Rpl13', 'Rps25', 'Rpl10-ps3', 'Rplp1', 'Rpl4', 'Rps27l',
                        'Rpl29', 'Rps27rt', 'Rpsa', 'Rpl14', 'Rps12', 'Rps15', 'Rpl41', 'Rps26',
                        'Rps27a', 'Rpl26', 'Rpl23a', 'Rps6kb1', 'Rpl23', 'Rpl19', 'Rpl27',
                        'Rpl38', 'Rps7', 'Rps29', 'Rpl36al', 'Rps6kl1', 'Rps6ka5', 'Rps23',
                        'Rpl15', 'Rps24', 'Rpl36a-ps1', 'Rpl37', 'Rpl30', 'Rpl8', 'Rpl3',
                        'Rps19bp1', 'Rpl39l', 'Rpl35a', 'Rpl24', 'Rps6ka2', 'Rps2', 'Rpl3l',
                        'Rps10', 'Rpl10a', 'Rps28', 'Rps18', 'Rpl7l1', 'Rpl36', 'Rpl36-ps4',
                        'Rps14', 'Rpl17', 'Rps6kb2', 'Rps6ka4', 'Rpl9-ps6', 'Rpl39', 'Rpl10',
                        'Rps4x', 'Rps6ka6', 'Rpl36a', 'Rps6ka3'),
           mt.genes = c('mt-Nd1', 'mt-Nd2', 'mt-Co1', 'mt-Co2', 'mt-Atp8', 'mt-Atp6', 'mt-Co3',
                        'mt-Nd3', 'mt-Nd4l', 'mt-Nd4', 'mt-Nd5', 'mt-Nd6', 'mt-Cytb'),
           hb.genes = c('Hbb-bt', 'Hbb-bs', 'Hbb-y', 'Hbs1l', 'Hba-a1', 'Hba-a2', 'Hbq1a',
                        'Hbegf'),
           sex.genes = c('Xist', 'Ddx3y', 'Eif2s3y', 'Erdr1', 'Gm29650', 'Kdm5d', 'Uty'),
           cc.genes = list('s.genes' =  c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6",
                                          "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2", "Rpa2", "Nasp", 
                                          "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2",
                                          "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2",
                                          "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8"),
                           'g2m.genes' = c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2",
                                           "Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2",
                                           "Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1",
                                           "Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2",
                                           "Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3",
                                           "Gas2l3","Cbx5","Cenpa" ))
         ))

config = new('config')
if(! dir.exists(config@result.path)){
  dir.create(config@result.path)
}
