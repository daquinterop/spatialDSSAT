# spatialDSSAT
Just a set of modules and functions to run DSSAT in spatial mode. With this package you will be able to simulate up to 99 different locations per run. You can define different weather, soil, nitrogen application, planting date, and cultivar for each location. Some simulation options can be defined. This allows the user to test different scenarios, like no water stress, no nitrogen stress, nitrogen fixation, etc.

This library is built on top of [DSSATTools](https://github.com/daquinterop/Py_DSSATTools). `run` Module contains the `GSRun` class. Same as `DSSAT` class hosts the single-treatment simulation environment in DSSATTools, GSRun hosts the multiple-location simulation environment in this package. Three steps are required to run an spatial simulation:

1. Create a GSRun object
2. Add treatments/locations
3. Run 

When GSRun is instantiated the crop name has to be passed as argument:

```python
gs = GSRun(crop_name="Maize")
```
After an instance of GSRun is created each location/treatment to be simulated is added using the `add_treatment` method:
```python
gs.add_treatment(
    soil="PATH/TO/.SOL",
    weather="PATH/TO/.WTH",
    nitro=[(0, 40), (30, 50)], # 40kg/ha at planting (0 dap), and 50kg/ha 30 dap
    planting=dateitme(2020, 1, 30),
    cultivar="990001"
)
```
Soil and weather are defined as the path to a `.WTH` and a `.SOL` file. Defining weather and soil as path to the DSSAT weather and soil files makes this tool independent of the soil or weather dataset used. This package also includes some tools to process netCDF files to create `.WTH` files. The `.SOL` file must contain only one soil profile. Nitrogen applications are defined as a list, with as many elements as nitrogen applications. Each application is defined by a tuple where the first element indicates the application date (days after planting) and the second element the nitrogen rate (kg/ha). It is assumed that the applied fertilizer is Urea. Planting date is defined as a datetime object, and cultivar is the cultivar code. Cultivar must be included in the DSSAT cultivar file.

After you've added the treatments/locations you want, you're set to run the model. The `run` method runs the model for the defined crop and treatments. It returns a pandas DataFrame with the results of the simulation.
```python
out = gs.run()
```
If you want to modify the simulation options you can pass them when running the model:
```python
out = gs.run(sim_controls={"WATER": "N"})
```
It the previous example the same set of treatments are run, but the processes that lead to water stress are not simulated. Similarly you can deactivate the processes that lead to nitrogen stress by setting `sim_controls={"NITRO": "N"}`.