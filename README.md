# BAC_firing

I N T R O D U C T I O N

Back-propagating action potential-activated Ca2+ spike (BAC) firing, first described by Larkum et al. (1999), is a burst of action potentials (APs) followed by a Ca2+-AP in the apical dendrite, resulting from the combined injection of current at the soma and apical dendrite. Due to backpropagation of APs, the threshold for a dendritic Ca2+ spike is subsequently lowered and leads to a burst firing. A combined stimulation of somatic input and apical EPSP-like input for a duration of 5ms will induce BAC firing (Larkum, 2013). 

This project replicates BAC firing using the reduced-morphology model by Mäki-Marttunen et al. (2018), and reviews the biophysical processes underlying BAC firing using some of the criteria presented by Bast and Oberlaender (2021) and Goetz et al. (2021). Furthermore, a Python script was written to detect such BAC criteria. This work aims to lay the groundwork for further developing a generalized BAC firing detection script that can be applied for other computational models.


I N S T R U C T I O N

Installing NEURON necessary for demonstrating simulation of any hoc file:
https://www.neuron.yale.edu/neuron/

Reduced-morphology model by Mäki-Marttunen et al. (2018) can be found here:
https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=187474#tabs-1

Marttunen_Template_CaFix.hoc serves merely as a template for developing simulations in hoc. It cannot be opened within NEURON interface. It describes the basic properties (channels and their features) and components (soma, dendrites, etc.) of a passive neuron.
Marttunen_Simulation.hoc is a coding script written in hoc to stimulate the neuron, i.e. turning it into an activen neuron which can be openend within NEURON. 



                                              (W I L L  B E  C O N T I N U E D)
                                              
                                              
                                              
R E F E R E N C E S

Bast, A., & Oberlaender, M. (2021). Ion channel distributions in cortical neurons are optimized for energy-efficient active dendritic computations. bioRxiv. https://doi.org/10.1101/2021.12.11.472235

Goetz, L., Roth, A., & Häusser, M. (2021). Active dendrites enable strong but sparse inputs to determine orientation selectivity. Proceedings of the National Academy of Sciences, 118(30), e2017339118. https://doi.org/10.1073/pnas.2017339118

Larkum, M., Zhu, J. & Sakmann, B. (1999). A new cellular mechanism for coupling inputs arriving at different cortical layers. Nature, 398(6725), 338–341. https://doi.org/10.1038/18686

Larkum, M. (2013). A cellular mechanism for cortical associations: an organizing principle for the cerebral cortex. Trends in Neurosciences, 36(3), 141-151. https://doi.org/10.1016/j.tins.2012.11.006

Mäki-Marttunen, T., Halnes, G., Devor, A., Metzner, C., Dale, A. M., Andreassen, O. A., Einevoll, G. T. (2018). A stepwise neuron model fitting procedure designed for recordings with high spatial resolution: Application to layer 5 pyramidal cells. Journal of Neuroscience Methods, 293, 264-283. https://doi.org/10.1016/j.jneumeth.2017.10.007
