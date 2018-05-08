############################
# Group 1
# Last edited: 3/27
# Eigenvector Centrality Values - Month 2
############################

library(igraph)

#Generate network using data from python script

graph = make_graph(~"Dennd2d"-"Igf2r",
"Kbtbd11"-"Arhgap9":"Slc2a4rg-ps":"Ccl5":"Cd28":"Gimap7":"Leprotl1":"Trav7d-3":"Gimap3":"cd28":"LOC385791":"LOC386513":"Trav9-1":"Trav7-1":"Tube1":"Xlr4c",
"Trav9d-3"-"Dgka":"Pik3r1",
"Rasal3"-"Prkch",
"Arhgap9"-"Kbtbd11":"Gimap7":"Trav7-1",
"Arl4c"-"Fxyd5":"Itgb2",
"Atp1b3"-"Bcl11b":"Slc2a4rg-ps":"Tes",
"Bcl11b"-"Atp1b3",
"Slc2a4rg-ps"-"Kbtbd11":"Atp1b3":"Ccl5":"Leprotl1":"cd28":"LOC385081":"Tes",
"C920016N10Rik"-"Cd163":"Slc40a1":"Vcam1",
"Ccl5"-"Kbtbd11":"Slc2a4rg-ps":"Cd28":"Fam105a":"Gimap7":"cd28":"LOC385081":"Trav9-1",
"Cd163"-"C920016N10Rik":"Satb1":"Slc40a1":"Thy1",
"Cd27"-"Ppic",
"Cd28"-"Kbtbd11":"Ccl5":"Gimap7":"Leprotl1":"cd28":"LOC386513":"Trav7d-3":"Trav9-1":"Trav7-1":"Sema4f",
"Cd3d"-"Cd3e":"Cd3g":"Cd6":"Fxyd5":"Lck":"Tcf7":"Trbv13-2",
"Cd3e"-"Cd3d":"Cd3g":"Dgka":"Lck":"Tcf7",
"Cd3g"-"Cd3d":"Cd3e":"Peli1":"Trbv13-2",
"Cd6"-"Cd3d":"Itgb7",
"Cd8b1"-"Selplg":"Slfn1",
"D12Ertd551e"-"Pik3r1",
"Dgka"-"Trav9d-3":"Cd3e":"Lef1":"Pik3r1",
"Dpysl2"-"Selplg",
"Faah"-"Leprotl1",
"Fam102a"-"Tes",
"Fam105a"-"Ccl5":"Sgk1":"Txk",
"Fxyd5"-"Arl4c":"Cd3d":"Itgb2":"Trbv13-2":"Tmsb10",
"Gimap4"-"LOC100041103":"Slc9a3r1",
"Gimap7"-"Kbtbd11":"Arhgap9":"Ccl5":"Cd28":"Trav7d-3":"cd28":"LOC385791":"LOC386513":"Trav9-1":"Trav7-1",
"Hcst"-"Ms4a6b",
"Igf2r"-"Dennd2d",
"Il18r1"-"LOC100046087":"LOC385791":"Rnf125",
"Itgb2"-"Arl4c":"Fxyd5",
"Itgb7"-"Cd6":"Trbv13-2",
"Itk"-"LOC100046087",
"Lat"-"Leprotl1":"LOC100046087":"cd28":"Rnf125":"Xlr4c",
"Lck"-"Cd3d":"Cd3e",
"Lcp2"-"Sgk1",
"Ldb1"-"LOC236170":"LOC385081":"scl000349.1_0",
"Lef1"-"Dgka":"scl000349.1_0",
"Leprotl1"-"Kbtbd11":"Slc2a4rg-ps":"Cd28":"Faah":"Lat":"cd28":"LOC386513":"Rnf125":"Tes":"Trib2":"Xlr4c",
"LOC100040243"-"Tmem66",
"LOC100041103"-"Gimap4",
"Trav7d-3"-"Kbtbd11":"Gimap7":"LOC385791":"LOC386513":"Trav7d-3":"Trav9-1":"Trav7-1":"Sema4f",
"Gimap3"-"Kbtbd11":"Ppm1b",
"LOC100046087"-"Il18r1":"Itk":"Lat":"Xlr4c",
"cd28"-"Kbtbd11":"Slc2a4rg-ps":"Ccl5":"Cd28":"Gimap7":"Lat":"Leprotl1":"LOC386513":"Trav7d-3":"Trav9-1":"Trav7-1":"Rnf125":"Xlr4c",
"LOC236170"-"Ldb1":"LOC385081":"scl000349.1_0",
"LOC385081"-"Slc2a4rg-ps":"Ccl5":"Ldb1":"LOC236170":"LOC386513":"Trav9-1":"scl000349.1_0",
"LOC385791"-"Kbtbd11":"Gimap7":"Il18r1":"Trav7d-3":"Trav9-1":"Tube1",
"LOC386513"-"Kbtbd11":"Cd28":"Gimap7":"Leprotl1":"Trav7d-3":"cd28":"LOC385081":"Trav7d-3":"Trav9-1":"Trav7-1":"scl000349.1_0",
"Trav7d-3"-"Cd28":"Trav7d-3":"cd28":"LOC386513":"Trav9-1":"Trav7-1":"Sema4f":"Trat1",
"Trav9-1"-"Kbtbd11":"Ccl5":"Cd28":"Gimap7":"Trav7d-3":"cd28":"LOC385081":"LOC385791":"LOC386513":"Trav7d-3":"Trav7-1":"scl000349.1_0":"Sema4f",
"Trav7-1"-"Kbtbd11":"Arhgap9":"Cd28":"Gimap7":"Trav7d-3":"cd28":"LOC386513":"Trav7d-3":"Trav9-1":"Trav6d-4":"Sema4f",
"Trav6d-4"-"Trav7-1",
"Ms4a6b"-"Hcst",
"Peli1"-"Cd3g",
"Pik3r1"-"Trav9d-3":"D12Ertd551e":"Dgka",
"Ppic"-"Cd27":"Tcf7",
"Ppm1b"-"Gimap3",
"Prkch"-"Rasal3",
"Rnf125"-"Il18r1":"Lat":"Leprotl1":"cd28":"Xlr4c",
"Satb1"-"Cd163":"Sepp1":"Slc40a1":"Thy1",
"scl000349.1_0"-"Ldb1":"Lef1":"LOC236170":"LOC385081":"LOC386513":"Trav9-1",
"Selplg"-"Cd8b1":"Dpysl2",
"Sema4f"-"Cd28":"Trav7d-3":"Trav7d-3":"Trav9-1":"Trav7-1":"Trat1",
"Sepp1"-"Satb1",
"Sgk1"-"Fam105a":"Lcp2",
"Slc40a1"-"C920016N10Rik":"Cd163":"Satb1":"Vcam1",
"Slc9a3r1"-"Gimap4",
"Slfn1"-"Cd8b1",
"Tcf7"-"Cd3d":"Cd3e":"Ppic",
"Trbv13-2"-"Cd3d":"Cd3g":"Fxyd5":"Itgb7",
"Tes"-"Atp1b3":"Slc2a4rg-ps":"Fam102a":"Leprotl1",
"Thy1"-"Cd163":"Satb1",
"Tmem66"-"LOC100040243",
"Tmsb10"-"Fxyd5",
"Tnfaip3"-"Traf1":"Xlr4c",
"Traf1"-"Tnfaip3":"Xlr4a",
"Trat1"-"Trav7d-3":"Sema4f",
"Trib2"-"Leprotl1",
"Tube1"-"Kbtbd11":"LOC385791",
"Txk"-"Fam105a",
"Vcam1"-"C920016N10Rik":"Slc40a1",
"Whsc1l1"-"Xlr4c",
"Xlr4a"-"Traf1",
"Xlr4c"-"Kbtbd11":"Lat":"Leprotl1":"LOC100046087":"cd28":"Rnf125":"Tnfaip3":"Whsc1l1")



# Visualize data
plot(graph)
# Get eigenvector centrality values
eigenvector_centrality_values <- eigen_centrality(graph)


#Create a column of the eigenvector centrality values
x <- matrix(eigenvector_centrality_values$vector, nrow=length(eigenvector_centrality_values$vector))
#Create a column of the names of the nodes
y <- matrix(names(eigenvector_centrality_values$vector), nrow=length(eigenvector_centrality_values$vector))


#combine to columns together
final <- cbind(x,y)


#open file 
file_output <- file("m2_evc_values.txt")

#write table to file
write.table(final, "m2_evc_values.txt", sep="\t")

close(file_output)


