import strformat
import strabo/lidar
import kdtree
import math
import cligen

##nim c -d:danger.\"program name"
##.\"program name".exe getinfo

#####Sithole algorithm based off a slope based filter proposed by Vosselman, (2000) and modified to include tophat morphological filter by Sithole, (2001).  The algorithm takes user input variables and filters any .las file based on slope.#####
#####Nathan Adams, ID 0928773#####
proc run(wd: string = "C:/Users/Natha/Data/", lidarFile: string = "1km_input.las", outputFile: string = "Demo.las", searchDist, minOffGroundPtHeight, maxDownwardSlope: float) = 

    let
        searchDistSqr = searchDist * searchDist #user defined 
    echo("Reading Points...")

    var
        las = newLasFile(lidarFile) 
        output = initializeUsingFile(outputFile, las) 
        points = newSeqOfCap[array[K, float]](las.header.numberOfPoints) 
        values = newSeqOfCap[int](las.header.numberOfPoints) 
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

    for p in 0..<las.header.numberOfPoints.int:
        pd = las.getPointRecord(p)
        maximumZ = minElevation[p]
        if pd.classField.classification != lowPoint and pd.classField.classification != highNoise:
            var ret = tree.withinRadius([pd.x, pd.y], radius=searchDistSqr, sortResults=false)
            for (_, n, dist) in ret:
                if minElevation[n] > maximumZ:
                    maximumZ = minElevation[n]      #If the minElevation is greater than the hgihest elevation in a neighbourhood
                                                    #maximumZ is increased to match that of minElevation
        tophatz.add(pd.z - maximumZ) 

        progress = (100.0 * (p.float + 1.0) / las.header.numberOfPoints.float).int  
        if progress mod 5 == 0 and progress != oldProgress:
            oldProgress = progress
            echo(&"Dilation: {progress}%")
        
    #Slope-based filter
    for p in 0..<las.header.numberOfPoints.int: 
        pd = las.getPointRecord(p)
        if pd.bitField.returnNumber() == pd.bitField.numberOfReturns():
            var ret = tree.withinRadius([pd.x, pd.y], radius=searchDistSqr, sortResults=false) 
            for (_, n, dist) in ret:                                                     
                run = sqrt(dist)
                pdNeighbour = las.getPointRecord(n)                                                                          
                if pdNeighbour.classField.classification != lowPoint and pdNeighbour.classField.classification != highNoise: 
                    rise = tophatz[p] - tophatz[n]
                    slope = radToDeg(arctan(rise/run))
                    if slope > maxDownwardSlope and rise >= minOffGroundPtHeight:
                        offGroundPoint[p] = true #if our point meets all of the above constraints it will be changed in the boolean array as a "true" value, signifying an offgroundpoint

        progress = (100.0 * (p.float + 1.0) / las.header.numberOfPoints.float).int 
        if progress mod 5 == 0 and progress != oldProgress:
            oldProgress = progress
            echo(&"Slope-based filter: {progress}%")
        
    #Output pointrecord
    for p in 0..<las.header.numberOfPoints.int: 
        if offGroundPoint[p] == false:  
            pd = las.getPointRecord(p)
            output.addPointRecord(pd)

    output.write()
    

#Beginning of user-defined variables using Cligen#
proc getinfo() = 
    echo("/////////////////////////////////")                 #C:\Users\Natha\Data\GEOG6060Data
    echo("//Sithole + White Tophat Filter//")                 #Guelph.las
    echo("/////////////////////////////////")                 #GuelphDemo().las

    write(stdout, "Working Directory: ") 
    let wd = readLine(stdin)

    write(stdout, "LiDAR Input File: ") 
    let lidarFile = readLine(stdin)
    doAssert len(lidarFile) > 0, "Error input file name is empty."

    write(stdout, "Output File: ") 
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

    #Appending the file name to the working directory if entered correctly#
    var inFile = if wd.len() > 0:
            wd & lidarFile 
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
main()

