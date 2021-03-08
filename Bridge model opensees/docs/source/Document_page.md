# Pybridge - Opensees Module

This page contains the guidelines 
for using the Opensees (OP) module for grillage
analysis.

   

## `Bridge` class

The ```Bridge``` class object contains all information from bridge pickle file
and hold it. 

The class then creates an ```OpenseesModel``` object which uses ```Openseespy``` methods to 
create the bridge model within the ```Openseespy``` framework.

The ```Bridge``` class object and ```OpenseesModel``` object is only instantiated through the ```Grillage``` class

The ```OpenseesModel``` object contains internal functions that communicates with ```Openseespy``` framework

Example: Using the bridge class
____________________

.. code-block::python

    # initialize Bridge class object within Grillage class instance
    self.OPBridge = OpenseesModel(self.bridgepickle["Nodedetail"], self.bridgepickle["Connectivitydetail"],
                                        self.bridgepickle["beamelement"], self.bridgepickle["Memberdetail"],
                                      self.bridgepickle["Member transformation"])
    # assign properties of concrete and steel
    self.OPBridge.assign_material_prop(self.bridgepickle["concreteprop"], self.bridgepickle["steelprop"])
    # send attribute to OP framework to create OP model
    self.OPBridge.create_Opensees_model()

    # time series and load pattern options
    self.OPBridge.time_series()
    self.OPBridge.loadpattern()


## ```Grillage``` class

The ```Grillage``` class takes two inputs: (1) a bridge pickle file, and (2) a ```vehicle``` class object.

Example: Using the grillage class
____________________

.. code-block::python

    # Properties of truck
    axlwts = [800, 3200, 3200]
    axlspc = [7, 7]
    axlwidth = 5
    initial_position = [0, 3.0875]
    travel_length = 50
    increment = 2
    direction = "X"

    # create Truck object
    RefTruck = vehicle(axlwts, axlspc, axlwidth, initial_position, travel_length, increment, direction)

    # Open bridge pickle file
    with open("save.p","rb") as f:
        refbridge = pickle.load(f)
    
    # create Grillage object
    RefBridge = Grillage(refbridge,RefTruck)
    
    # perform moving truck analysis
    RefBridge.perfromtruckanalysis()

## Guidelines for bridge pickle file

A bridge model is loaded through a pickle
file which contain the bridge information.

The pickle file is formatted to communicate with the classes and methods of the module. Note: no modifications should be 
 performed with respect to the bridge pickle fill unless reviewed with changes in the classes and methods. 
 
The pickle file contains a dictionary of dataframas (pandas), each details the properties
of the bridge. The entries include:
1) Node data
2) Connectivity data
3) Member properties
4) Material properties

For more information of the data frame and its content please refer to the 
example excel file - ReferenceBridge.xlsx

Example: bridge pickle file
____________________

