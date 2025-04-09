Run the build.sh file to build the libraries and the application.
The run ./build/fv99Hpc 
Then copy the file ./build/segment.vtp and try opening it with paraview to see if it's worked. It should be a line segment.

# Run examples

To get help with the application

```
./fv99 -h
```


To run on a small toy dataset
```
./fv99 -f ../data/data.vtu --outputSheetPolygons ./toyExampleSheetPolygons.txt --outputSheetFibersFolder ./toyExampleFibers --fiberSampling 20 --sheetOutputCount 30 -e 0.001
```

Explanation of the parameter

Text file with the sheets of the Reeb space as polygons in plain text
```
--outputSheetPolygons ./toyExampleSheetPolygons.txt
```

Representative fibers for the largest sheets in terms of area
```
--outputSheetFibersFolder ./toyExampleFibers
```

How many fibers to sample per sheet
```
--fiberSampling 20
```

How many sheets to sample
```
--sheetOutputCount 30
```


To run on a slighly bigger data set

```
./fv99 -f ../data/downsample-15-875.vtu --outputSheetPolygons ./downsample15.txt --outputSheetFibersFolder ./downsample15Fibers --fiberSampling 20 --sheetOutputCount 30 -e 0.001

```

Big boy
```
./fv99 -f ~/Projects/data/reeb-space-test-data/ttk/downsample-4-46800.vtu --outputSheetFibersFolder ./ttk-downsample-4 --fiberSampling 20 --sheetOutputCount 30 -e 1e-2

```
