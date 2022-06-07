// Fiji Macro to make kymographs
// First, install kymograph plugin from http://biop.epfl.ch/TOOL_KYMOGRAPH.html which was
// developed in the paper Processive movement of single kinesins on crowded
// microtubules visualized using quantum dots. Seitz A, Surrey T., EMBO J. 2006
// Jan 25;25(2):267-77. Epub 2006 Jan 12.
for(i=1;i<31;i++){
open("C:/XX/simulation_result/prom2_"+i+".tiff");
open("C:/XX/simulation_result/prom2.roi");
run("KymographMax", "linewidth=11");
saveAs("Tiff", "C:/XX/simulation_result/kymo_prom2_"+i+"_noise.tiff");
close();
close();
}
