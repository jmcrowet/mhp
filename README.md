# mhp
Pymol plugin for Mean Hydrophobicity Potential calculations

## DESCRIPTION

  This plugin computes the MHP (Mean Hydrophobicity Potential) values 
  for the specified pdb file and draw the isopotentials in PyMol.
  MHP values are computed for each positions on a user defined grid and 
  saved in a *.dx file. 
  MHP isopotentials are distant from the molecular surface so you have to 
  define a distance between the molecule and the box walls.
  
## USAGE 

  Through the interface

     If you have not a *.dx file yet, use the MHP Compute procedure 
     of this plugin. Fill in the parameters (default values are 
     generally OK) and click on OK.
     MHP computaions take some time so if you have *.dx file, 
     you should use the MHP Display procedure.

  On the command line

     MHP_pymol *.pdb, *.dx, [d, delta, iso_phi, iso_pho]

        Terms in brackets are optional and their default values 
        are d = 5, delta = 1, iso_phi = 0.1, iso_pho = -0.1
        d correspond to the distance between the solute and the box (A)
        delat correspond to the mesh size (A)
        iso_phi and iso_pho correspond to the isopotential values

     MHP_display *.dx [iso_phi, iso_pho]

        MHP_display draw the isopotentials in PyMol from the *.dx file.

     MHP_Help

Created by Jean-Marc Crowet (jmcrowet@ulg.ac.be)
CBMN - University of Liege (www.fsagx.ac.be/bp/)

