# D<sup>2</sup><sub>min</sub>(t, $\Delta$ t)  
computed D<sup>2</sup><sub>min</sub> by Langer & Falk, see:  
* Falk, M. L., & Langer, J. S. (1998). Dynamics of viscoplastic deformation in amorphous solids. Physical Review E, 57(6), 7192-7205. doi:DOI 10.1103/PhysRevE.57.7192  
* Cubuk, E. D., Ivancic, R. J. S., Schoenholz, S. S., Strickland, D. J., Basu, A., Davidson, Z. S., . . . Liu, A. J. (2017). Structure-property relationships from universal signatures of plasticity in disordered solids. Science, 358(6366), 1033-1037. doi:doi:10.1126/science.aai8830
  
<boost/filesystem.hpp> in boost library is used to create file & dirs,  
  
Eigen3 library is used to carry out matrix computation  

files:
* **configurtaion.h** contains the basic *Configuration* class  
* **D2.h** contains *Configuration_neighbours* class and *D2* class  
     * *Configuration_neighbours* is derived from *Configuration*, which contains neighbour lists of particles  
    * *D2* is used to compute D<sup>2</sup><sub>min</sub>(t, $\Delta$ t) 
* **input.h**, class *Input*, a class used to read files as line by line
* **particle.h**, contains sturcture *Particle* 

* dir **PythonTools** contains some automation tools and useful tools like *isoconfiguration average*(see: Widmer-Cooper, A., Harrowell, P., & Fynewever, H. (2004). How Reproducible Are Dynamic Heterogeneities in a Supercooled Liquid? Physical Review Letters, 93(13), 135701. doi:10.1103/PhysRevLett.93.135701
)

