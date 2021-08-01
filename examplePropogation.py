# examplePropogation.py



from loadTLE import *
from satFunctions import *
from locations import *

from skyfield import api


#Important things
filename = "45748.txt"
location = locations["Lemmon"]
timeStart = dt.datetime(2020,9,8,3,44,47,760000, tzinfo=api.utc)
timeEnd = dt.datetime(2020,9,8,3,44,47+3,760000, tzinfo=api.utc)


#Load TLE from a local file
#Normally loads a whole list fo TLEs so grab the first one for now
tle = loadTLEFile(filename)[0]
printTLE(tle)

#Compute the start points and end points the the streak
start = computeEphemeris(tle, location, timeStart)
end = computeEphemeris(tle, location, timeEnd)


#Print results
printEphemeris(start)
printEphemeris(end)





