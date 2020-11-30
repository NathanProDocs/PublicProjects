import strformat
import strabo/lidar
import kdtree
import math
import cligen
import strutils

proc run(wd: string = "C:/Users/Natha/Data/", lidarFile: string = "1km_input.las", outputFile: string = "Demo.las", searchDist, minOffGroundPtHeight, maxDownwardSlope: float) = 

    let
        searchDistSqr = searchDist * searchDist #user defined 
    echo("Reading Points...")

    var
        las = newLasFile(lidarFile) #Reading in the initial LiDAR point cloud
        output = initializeUsingFile(outputFile, las) #Creates the new las file
        points = newSeqOfCap[array[K, float]](las.header.numberOfPoints) #creates an array of points
        values = newSeqOfCap[int](las.header.numberOfPoints) #creates an array of values for those points
        pd, pdNeighbour: PointRecord3 #gets the point record data from each point
        progress: int 
        old_progress = 1 
        offGroundPoint = newSeqOfCap[bool](las.header.numberOfPoints) 
        run: float
        slope: float
        rise: float

    for p in 0..<las.header.numberOfPoints.int:
        pd = las.getPointRecord(p)  
        if pd.bitField.returnNumber() == pd.bitField.numberOfReturns():   # check if the return number equals the total number of returns (e.g. if it is the 5th return, was there 5 returns total)
            points.add([pd.x, pd.y])                         
            values.add(p)  
            
            offGroundPoint.add(false)                                   
        else:
            offGroundPoint.add(true)                                      

        progress = (100.0 * (p.float + 1.0) / las.header.numberOfPoints.float).int  
        if progress mod 5 == 0 and progress != oldProgress:
            oldProgress = progress
            echo(&"Reading input data: {progress}%")

    echo("Building Kd-tree...")
    var tree = newkdtree[int](points, values)   

    echo("Finished building tree...") 

    #Below is white tophat transform
    #Erosion
    var
        minElevation = newSeqOfCap[float](las.header.numberOfPoints)
        minimumZ: float

    for p in 0..<las.header.numberOfPoints.int:
        pd = las.getPointRecord(p)
        minimumZ = pd.z
        if pd.classfield.classification != lowPoint and pd.classField.classification != highNoise:
            var ret = tree.withinRadius([pd.x, pd.y], radius=searchDistSqr, sortResults = false)
            for (_, n, dist) in ret:
                if las[n].z < minimumZ and las.getPointRecord(n).classField.classification != lowPoint:
                    minimumZ = las[n].z #finds the lowest elevation point in the neighbourhood
        minElevation.add(minimumZ)

        progress = (100.0 * (p.float + 1.0) / las.header.numberOfPoints.float).int  
        if progress mod 5 == 0 and progress != oldProgress:
            oldProgress = progress
            echo(&"Erosion: {progress}%")

    #Dilation
    var
        tophatz = newseqofcap[float](las.header.numberOfPoints)
        maximumZ: float
        #nnDist = newseqofcap[float](las.header.numberOfPoints)

    for p in 0..<las.header.numberOfPoints.int:
        pd = las.getPointRecord(p)
        maximumZ = minElevation[p]
        if pd.classField.classification != lowPoint and pd.classField.classification != highNoise:
            var ret = tree.withinRadius([pd.x, pd.y], radius=searchDistSqr, sortResults=false)
            for (_, n, dist) in ret:
                if minElevation[n] > maximumZ:
                    maximumZ = minElevation[n] #If the minElevation is greater than the hgihest elevation in a neighbourhood
                                                    #maximumZ is increased to match that of minElevation
        tophatz.add(pd.z - maximumZ) #adds values to an array representing the point of interest - highest elevation point in the neighbourhood

        progress = (100.0 * (p.float + 1.0) / las.header.numberOfPoints.float).int  
        if progress mod 5 == 0 and progress != oldProgress:
            oldProgress = progress
            echo(&"Dilation: {progress}%")
        
    #3.	Visit each point in the LAS file, if the current value of the point index in 'offGroundPoint' is false then perform a withinRadius search at the site of the point's x,y coordinates. 
    for p in 0..<las.header.numberOfPoints.int: 
        pd = las.getPointRecord(p)
        if pd.bitField.returnNumber() == pd.bitField.numberOfReturns():
            var ret = tree.withinRadius([pd.x, pd.y], radius=searchDistSqr, sortResults=false) 
            for (_, n, dist) in ret:                                                     
                run = sqrt(dist)
                pdNeighbour = las.getPointRecord(n)                                                                          
                if pdNeighbour.classField.classification != lowPoint and pdNeighbour.classField.classification != highNoise: #loop through each point and discard those which are classed as lowpoints or with highnoise
                    #4.	For each neighbouring point, measure the slope between the neighbour and the point of interest at the center of the search. 
                    #The rise = difference in elevation between the neighbour and the point of interest, the run is the distance separating them. slope = arctan(rise/run)--remember to convert to degrees.
                    rise = tophatz[p] - tophatz[n]
                    slope = radToDeg(arctan(rise/run))
                    #5.	Find the maximum DOWNWARD slope (i.e. if you measure rise as PointOfInterest elevation - neighbour elevation then you are looking for the maximum positive slope value...otherwise you want the maximum negative value). 
                    #We essentially only want to know the slopes between the point of interest and neighbours that are lower than it.
                    if slope > maxDownwardSlope and rise >= minOffGroundPtHeight:
                        #6.	If the maximum downward slope is above a user-specified threshold value (in degrees) then the point of interest is an off-ground point and you should mark it as true in the 'offGroundPoint' array.
                        offGroundPoint.add(true) #if our point meets all of the above constraints it will be added back to our boolean array as a "true" value, signifying an offgroundpoint
        progress = (100.0 * (p.float + 1.0) / las.header.numberOfPoints.float).int 
        if progress mod 5 == 0 and progress != oldProgress:
            oldProgress = progress
            echo(&"Slope-based filter: {progress}%")

    for p in 0..<las.header.numberOfPoints.int: 
        if offGroundPoint[p] == false:  #sq brackets for index, round for function
            pd = las.getPointRecord(p)
            output.addPointRecord(pd)

    output.write()

proc getinfo() = 
    echo("/////////////////////////////////")
    echo("/Sithole + White Tophat Filter/")
    echo("/////////////////////////////////")

    write(stdout, "Working Directory: ") #   C:\Users\Natha\Data\
    let wd = readLine(stdin)

    write(stdout, "LiDAR Input File: ") #   1km_input.las
    let lidarFile = readLine(stdin)
    doAssert len(lidarFile) > 0, "Error input file name is empty."

    write(stdout, "Output File: ") #   Demo.las
    let outputFile = readLine(stdin)
    doAssert len(outputFile) > 0, "Error output file name is empty."

    var searchDist = float(1.5) 
    write(stdout, "Search distance (default = 1.5): ")
    let searchDistUD = readline(stdin)
    if len(searchDistUD) > 0:
        searchDist = parseFloat(searchDistUD)
    
    var minOffGroundPtHeight = float(0.5)
    write(stdout, "Minimum off-ground Point Height (default = 0.5): ")
    let minOffGroundPtHeightUD = readline(stdin)
    if len(minOffGroundPtHeightUD) > 0:
        minOffGroundPtHeight = parsefloat(minOffGroundPtHeightUD)

    var maxDownwardSlope = float(50.0)
    write(stdout, "Maximum Downward Slope Angle (default = 50.0): ")
    let maxDownwardSlopeUD = readline(stdin)
    if len(maxDownwardSlopeUD) > 0:
        maxDownwardSlope = parsefloat(maxDownwardSlopeUD)

    var inFile = if wd.len() > 0:
            wd & lidarFile #This part is required to join the input files and output files to the working directory, will throw an error if not added
        else:
            lidarFile

    var outFile = if wd.len() > 0:
            wd & outputFile
        else:
            outputFile

    run(wd, inFile, outFile, searchDist, minOffGroundPtHeight, maxDownwardSlope)

proc main() = 
    dispatchMulti([run, help={
            "wd": "Working Directory.",
            "lidarFile": "Name of your input .las file.",
            "outputFile": "Name of output .las file.",
            "searchDist": "Radius distance used in neighbourhood calculation.",
            "minOffGroundPtHeight": "Set the minimum off ground point height threshold.",
            "maxDownwardSlope": "Set the maximum downward slope threshold in degrees."}
        ],
        [getinfo]
        )
echo("Done!")
main()
