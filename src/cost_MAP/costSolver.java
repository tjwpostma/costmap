/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cost_MAP;

import org.w3c.dom.ls.LSOutput;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Dictionary;
import java.util.Hashtable;
import java.util.List;

public class costSolver {
//
    int[] rows = {0,-1,0,1,0,-1,-1,1,1};
    int[] cols = {0, 0,1,0,-1,-1,1,1,-1};

    /**
     * Parses ESRI ASCII raster header to extract grid metadata.
     * Reads the first 6 lines of an ASCII raster file and extracts dimensions,
     * geographic coordinates, cell size, and no-data value.
     *
     * @param path File path to the .asc raster file
     * @return Dictionary containing keys: {Columns, Rows, xllCorner, yllCorner, CellSize, NoData}
     *         Example: {NoData=-9999, CellSize=0.008333333333, xllCorner=-88.683333333822,
     *         Columns=220, Rows=174, yllCorner=37.641666667611}
     * @throws FileNotFoundException if raster file not found
     * @throws IOException if error reading file
     */
    public Dictionary getHeader(String path) throws FileNotFoundException, IOException {

//returns the info as a header dictionary {NoData=-9999, CellSize=0.008333333333, xllCorner=-88.683333333822, Columns=220, Rows=174, yllCorner=37.641666667611}
//        ArrayList header = new ArrayList();
        Dictionary header = new Hashtable();
        BufferedReader  br = new BufferedReader(new FileReader(path));
        //initiate reader
        String line = br.readLine();
        //define delimeters
        String delims = "[ ]+";
        //parse first line to get number of columns
        String[] columnLine = line.split(delims);
        int columns = Integer.parseInt(columnLine[1]);
        header.put("Columns", columns);
        //read next line and then get number of rows
        line = br.readLine();
        String[] rowLine = line.split(delims);
        int rows = Integer.parseInt(rowLine[1]);
        header.put("Rows", rows);
        //create empty matrix
        line = br.readLine();
        String[] xllCornerLine = line.split(delims);
        double xllCorner = Double.parseDouble(xllCornerLine[1]);
        header.put("xllCorner", xllCorner );
        line = br.readLine();
        String[] yllCornerLine = line.split(delims);
        double yllCorner = Double.parseDouble(yllCornerLine[1]);
        header.put("yllCorner", yllCorner );
        line = br.readLine();
        String[] cellSizeLine = line.split(delims);
        double cellSize = Double.parseDouble(cellSizeLine[1]);
        header.put("CellSize", cellSize);
        line = br.readLine();
        String[] noDataline = line.split(delims);
        int noData = Integer.parseInt(noDataline[1]);
        header.put("NoData", noData);

        return header;
    }

    /**
     * Loads raster grid data from ESRI ASCII file into a 2D matrix.
     * Reads data values after the 6-line header and populates a double[rows][cols] array.
     * WARNING: Large rasters (>10K x 10K) consume significant heap memory.
     *
     * @param headerInfo Header Dictionary from getHeader() containing Rows and Columns
     * @param path File path to the .asc raster file
     * @return 2D double array [rows][columns] containing grid values; no-data cells are -9999
     * @throws FileNotFoundException if raster file not found
     * @throws IOException if error reading file
     */
    public double[][] getDetails(Dictionary headerInfo, String path) throws FileNotFoundException, IOException {

        //returns the big number chunk below the header
        BufferedReader br = new BufferedReader(new FileReader(path));

        double [][] aoiMatrix = new double[(int) headerInfo.get("Rows")][(int) headerInfo.get("Columns")];
        //initiate reader, read lines in sequence
        String line = br.readLine();  //ncol
        //read next line and then get number of rows
        line = br.readLine();  //nrows

        line = br.readLine();  //xll corner

        line = br.readLine();  //yll corner

        line = br.readLine();  //cell size

        line = br.readLine();  //no data

        line = br.readLine();
        while (line != null) {
            for (double[] aoiMatrix1 : aoiMatrix) {
                String[] values = line.split("[ ]+");
                for (int j = 0; j < values.length; j++) {
                    aoiMatrix1[j] = Double.parseDouble(values[j]);

                }
                line = br.readLine();
            }
        }
        br.close();

        return aoiMatrix;

    }

    /**
     * Applies landcover classification weights to create baseline cost surface.
     * Reads NLCD landcover codes and maps them to weights from Datasets/weights/landcover.txt.
     * Optionally incorporates population density for urban/developed land codes (21-24).
     * Generates 9-cell kernels for each grid cell storing costs to traverse to neighbors.
     *
     * Landcover codes:
     * - 11=water, 12=evergreen forest, 21-24=developed/urban, 31=barren, 41-43=forest,
     *   52=shrub/scrub, 71=grassland, 81-82=pasture/cropland, 90/95=wetland
     *
     * Population density tiers (for codes 21-24):
     * - 0-5 pers/km²: weight 0.75 | 5-25: 1.5 | 25-100: 2.5 | 100+: 5.0
     * - If population not selected: all developed land weight 5.0
     *
     * @param isSelectedPop If true, incorporates population.asc density overlay for urban cells
     * @param path Path to landcover.asc raster
     * @return Dictionary<Integer cellId, double[9] costs> where cellId starts at 1 (not 0)
     * @throws Exception if weights file or landcover raster cannot be read
     */
    public Dictionary landcoverInput(boolean isSelectedPop, String path) throws FileNotFoundException, IOException, Exception {

        //Read in landcover weighting
        ArrayList weights = new ArrayList();
        BufferedReader br1 = new BufferedReader(new FileReader("Datasets/weights/landcover.txt"));
        String line = br1.readLine();
        while ((line = br1.readLine()) != null) {
            String[] splited = line.split("\\s+");
            weights.add(splited[2]);
        }
        br1.close();

        Dictionary headerInfo = getHeader(path);
        double[][] nlcdMatrix = getDetails(headerInfo, path);


        //Create Output matrix
        double[][] tempMatrix = new double[nlcdMatrix.length][nlcdMatrix[0].length];
        double[][] popMatrix = null;

        //Read in population file if it exists
        if (isSelectedPop == true) {
            System.out.println("Importing Population Data ...");
            popMatrix = getDetails(headerInfo, "Datasets/ASCII/population.asc");
        }

        // Precompute parameters once (outside the loops)
        int noData = (int) headerInfo.get("NoData");

        double w0 = Double.parseDouble((String) weights.get(0));
        double w1 = Double.parseDouble((String) weights.get(1));
        double w2 = Double.parseDouble((String) weights.get(2));
        double w3 = Double.parseDouble((String) weights.get(3));
        double w4 = Double.parseDouble((String) weights.get(4));
        double w5 = Double.parseDouble((String) weights.get(5));
        double w6 = Double.parseDouble((String) weights.get(6));
        double w7 = Double.parseDouble((String) weights.get(7));
        double w8 = Double.parseDouble((String) weights.get(8));

        for (int i = 0; i < nlcdMatrix.length; i++) {
            for (int j = 0; j < nlcdMatrix[0].length; j++) {
                int value = (int) nlcdMatrix[i][j];
                switch(value){
                    case 11:
                        tempMatrix[i][j] = w0;
                        break;
                    case 12:
                        tempMatrix[i][j] = w1;
                        break;
                    case 21:
                    case 22:
                    case 23:
                    case 24:
                        if (isSelectedPop) {
                            double cellPop = cellPop(headerInfo, popMatrix, i, j);

                            if (cellPop == 0) {
                                tempMatrix[i][j] = 0.75;
                            } else if (cellPop <= 5) {
                                tempMatrix[i][j] = 1.0;
                            } else if (cellPop <= 25) {
                                tempMatrix[i][j] = 1.5;
                            } else if (cellPop <= 100) {
                                tempMatrix[i][j] = 2.5;
                            } else {
                                tempMatrix[i][j] = 5.0;
                            }
                        } else {
                            tempMatrix[i][j] = 5.0;
                        }
                        break;
                    case 31:
                        tempMatrix[i][j] = w2;
                        break;
                    case 41:
                    case 42:
                    case 43:
                        tempMatrix[i][j] = w3;
                        break;

                    case 52:
                        tempMatrix[i][j] = w4;
                        break;

                    case 71:
                        tempMatrix[i][j] = w5;
                        break;
                    case 81:
                        tempMatrix[i][j] = w6;
                        break;
                    case 82:
                        tempMatrix[i][j] = w7;
                        break;
                    case 90:
                    case 95:
                        tempMatrix[i][j] = w8;
                        break;
                    default:
                        tempMatrix[i][j] = noData;

                }

        }
        }
        Dictionary costList = new Hashtable();
        int cell = 1;  // Cells have 1-based indexing in cost list
        for (int i = 0; i < tempMatrix.length; i++) {
            for (int j = 0; j < tempMatrix[0].length; j++) {
                    double[] landKernel = kernel(tempMatrix, i, j);
                    double[] costs = solveLand(landKernel);
                    costList.put(cell, costs);
                    cell+=1;
            }
        }


       return costList;
    }

    /**
     * Alternative landcover processing that applies federal land overlay AFTER landcover weighting.
     * Multiplies landcover costs by federal land agency multipliers from fed.asc overlay.
     * Use this mode when federal land restrictions should override standard landcover weights.
     *
     * Federal agency codes and multipliers in fed.asc:
     * - 1 (BLM), 2 (BOR): 0.5x | 3 (DOD): 2x | 4 (FS): 1.5x | 5 (FWS): 2.5x
     * - 6 (NPS): 2x | 7 (Other): 3x | 8 (TVA): 0.75x | 9 (State Parks): 2x | 10 (Reservations): 50x
     *
     * Currently only one mode (landcoverInput or landRowInput) is used per run;
     * GUI checkbox would need to be added to toggle between modes.
     *
     * @param isSelectedPop If true, incorporates population.asc density overlay for urban cells
     * @param path Path to landcover.asc raster
     * @return Dictionary<Integer cellId, double[9] costs> with federal land multipliers applied
     * @throws Exception if weights file or rasters cannot be read
     */
    public Dictionary landRowInput(boolean isSelectedPop, String path) throws FileNotFoundException, IOException, Exception {

        //Read in landcover weighting
        ArrayList weights = new ArrayList();
        BufferedReader br1 = new BufferedReader(new FileReader("Datasets/weights/landrows.txt"));
        String line = br1.readLine();
        while ((line = br1.readLine()) != null) {
            String[] splited = line.split("\\s+");
            weights.add(splited[2]);
        }
        br1.close();

        Dictionary headerInfo = getHeader(path);
        double[][] nlcdMatrix = getDetails(headerInfo, path);
        double[][] fedMatrix = getDetails(headerInfo, "Datasets/ASCII/fed.asc");


        //Create Output matrix
        double[][] tempMatrix = new double[nlcdMatrix.length][nlcdMatrix[0].length];
        double[][] popMatrix = null;

        //Read in population file if it exists
        if (isSelectedPop == true) {
            System.out.println("Importing Population Data ...");
            popMatrix = getDetails(headerInfo, "Datasets/ASCII/population.asc");
        }

        // Precompute parameters once (outside the loops)
        int noData = (int) headerInfo.get("NoData");

        double w0 = Double.parseDouble((String) weights.get(0));
        double w1 = Double.parseDouble((String) weights.get(1));
        double w2 = Double.parseDouble((String) weights.get(2));
        double w3 = Double.parseDouble((String) weights.get(3));
        double w4 = Double.parseDouble((String) weights.get(4));
        double w5 = Double.parseDouble((String) weights.get(5));
        double w6 = Double.parseDouble((String) weights.get(6));
        double w7 = Double.parseDouble((String) weights.get(7));
        double w8 = Double.parseDouble((String) weights.get(8));

        for (int i = 0; i < nlcdMatrix.length; i++) {

            for (int j = 0; j < nlcdMatrix[0].length; j++) {

                int value = (int) nlcdMatrix[i][j];
                switch(value){
                    case 11:
                        tempMatrix[i][j] = w0;
                        break;
                    case 12:
                        tempMatrix[i][j] = w1;
                        break;
                    case 21:
                    case 22:
                    case 23:
                    case 24:
                        if (isSelectedPop) {
                            double cellPop = cellPop(headerInfo, popMatrix, i, j);

                            if (cellPop == 0) {
                                tempMatrix[i][j] = 0.75;
                            } else if (cellPop <= 5) {
                                tempMatrix[i][j] = 1.0;
                            } else if (cellPop <= 25) {
                                tempMatrix[i][j] = 1.5;
                            } else if (cellPop <= 100) {
                                tempMatrix[i][j] = 2.5;
                            } else {
                                tempMatrix[i][j] = 5.0;
                            }
                        } else {
                            tempMatrix[i][j] = 5.0;
                        }
                        break;
                    case 31:
                        tempMatrix[i][j] = w2;
                        break;
                    case 41:
                    case 42:
                    case 43:
                        tempMatrix[i][j] = w3;
                        break;

                    case 52:
                        tempMatrix[i][j] = w4;
                        break;

                    case 71:
                        tempMatrix[i][j] = w5;
                        break;
                    case 81:
                        tempMatrix[i][j] = w6;
                        break;
                    case 82:
                        tempMatrix[i][j] = w7;
                        break;
                    case 90:
                    case 95:
                        tempMatrix[i][j] = w8;
                        break;
                    default:
                        tempMatrix[i][j] = noData;

                }

            }
        }

        for (int i = 0; i < fedMatrix.length; i++) {

            for (int j = 0; j < fedMatrix[0].length; j++) {

                int fed = (int) fedMatrix[i][j];
                double value = tempMatrix[i][j];

                switch (fed) {
                    case 1: // BLM
                    case 2: // BOR
                        tempMatrix[i][j] = value * 0.5;
                        break;
                    case 3: // DOD
                        tempMatrix[i][j] = value * 2.0;
                        break;
                    case 4: // FS
                        tempMatrix[i][j] = value * 1.5;
                        break;
                    case 5: // FWS
                        tempMatrix[i][j] = value * 2.5;
                        break;
                    case 6: // NPS
                        tempMatrix[i][j] = value * 2.0;
                        break;
                    case 7: // Other
                        tempMatrix[i][j] = value * 3.0;
                        break;
                    case 8: // TVA
                        tempMatrix[i][j] = value * 0.75;
                        break;
                    case 9: // State Parks
                        tempMatrix[i][j] = value * 2.0;
                        break;
                    case 10: // Reservations
                        tempMatrix[i][j] = value * 50.0;
                        break;
                    default:
                        tempMatrix[i][j] = noData;
                }
            }
        }

        Dictionary costList = new Hashtable();
        int cell = 1;  // Cells have 1-based indexing in cost list
        for (int i = 0; i < tempMatrix.length; i++) {
            for (int j = 0; j < tempMatrix[0].length; j++) {
                    double[] landKernel = kernel(tempMatrix, i, j);
                    double[] costs = solveLand(landKernel);
                    costList.put(cell, costs);
                    cell+=1;
            }
        }
       return costList;
    }

    /**
     * Converts raster population count to population density (persons/km²).
     * Accounts for latitude effects on longitude distance using haversine formula.
     * Used to apply density-based modifiers to urban landcover cells.
     *
     * CRITICAL: Verify lat/lon corner order matches raster orientation (yllCorner typically bottom-left).
     * Inverted corners will produce inverted/incorrect density calculations.
     *
     * @param headerInfo Header Dictionary from getHeader() containing yllCorner, xllCorner, CellSize
     * @param array Population raster matrix loaded via getDetails()
     * @param i Row index of cell
     * @param j Column index of cell
     * @return Population density in persons/km²
     */
    private double cellPop(Dictionary headerInfo,double[][] array, int i, int j) {
        double pop = array[i][j];
        double lon = (Double) headerInfo.get("xllCorner");
        double lat = (Double) headerInfo.get("yllCorner") + (((int) headerInfo.get("Rows") - (i+1)) * (Double) headerInfo.get("CellSize")) + ((Double) headerInfo.get("CellSize")/2);
        double cHeight = haversineDistance(lat, lon, lat + (Double) headerInfo.get("CellSize"), lon);
        double cWidth = haversineDistance(lat, lon, lat, lon + (Double) headerInfo.get("CellSize"));
        double cArea = cWidth * cHeight;
        return pop / cArea;
    }

    private double haversineDistance(double lat1, double lon1, double lat2, double lon2) {

        double earthRadius = 6369.15;
        double latDistance = toRad(lat2 - lat1);
        double lonDistance = toRad(lon2 - lon1);
        double a = (Math.sin(latDistance / 2) * Math.sin(latDistance / 2)
                + Math.cos(toRad(lat1)) * Math.cos(toRad(lat2))
                * Math.sin(lonDistance / 2) * Math.sin(lonDistance / 2));
        double c = (2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a)));
        double distance = earthRadius * c;
        return distance;
    }

    //convert to radians
    private double toRad(double latChange) {
        return (double) (latChange * Math.PI / 180);
    }


    /**
     * Applies geographic distance multipliers to cost network.
     * Scales costs based on actual cell-to-cell distances (not just cell count).
     * Different multipliers for cardinal (E/W/N/S) vs diagonal directions.
     * Must be called after landcoverInput() to have valid cost Dictionary.
     *
     * @param headerInfo Header Dictionary from getHeader()
     * @param distMult Distance multiplier matrix from distanceMultiplier() [rows][4]
     *                 where [i][0]=E/W, [i][1]=N/S, [i][2]=NE/SW, [i][3]=NW/SE
     * @param costList Dictionary from landcoverInput() with initial costs
     * @param path Path to landcover.asc (used to iterate through grid dimensions)
     * @return Updated costList with distance-scaled values
     */
    public Dictionary solveDistance(Dictionary headerInfo, double[][] distMult, Dictionary costList, String path) throws IOException {
        double[][] nlcdMatrix = getDetails(headerInfo, path);

        int indexNew = 1;
        for (int z = 0; z < nlcdMatrix.length; z++) {
            for (int j = 0; j < nlcdMatrix[0].length; j++) {
                if (indexNew < costList.size()) {
                    double[] costs = (double[]) costList.get(indexNew);
                    for (int i = 0; i < costs.length; i++) {
                        if (i == 2 || i == 4) {
                            costs[i] = round((costs[i]) * distMult[z][0], 2);
                        } else if (i == 1 || i == 3) {
                            costs[i] = round((costs[i]) * distMult[z][1], 2);
                        } else if (i == 5 || i == 6) {
                            costs[i] = round((costs[i]) * distMult[z][2], 2);
                        } else {
                            costs[i] = round((costs[i]) * distMult[z][3], 2);
                        }

                        costList.put(indexNew, costs);
                    }

                }
                indexNew+=1;
            }
        }
        return costList;
    }


    /**
     * Adds terrain cost increments based on slope angle and aspect direction.
     * Slope increments: 0.1 (0.5-1°), 0.2 (1-1.5°), 0.3 (1.5-2°), 0.4 (2-2.5°), 0.5 (2.5-3°), 1.0 (>3°).
     * Aspect-aware adjustments: +0.2 if upslope in direction of travel, +0.1 if crossslope.
     * Average new slope cost with center and neighbor per cell: (slopeCenter + slopeNeighbor) / 2.
     *
     * Aspect codes (from aspectInput()): 1-8 for N/NW/W/SW/S/SE/E/NE.
     * Kernel indexing for directions: 1=N, 2=E, 3=S, 4=W, 5=NW, 6=NE, 7=SE, 8=SW.
     *
     * @param costList Dictionary from landcoverInput() with baseline costs
     * @param isSelectedAspect If true, applies directional upslope penalties via aspect.asc
     * @param path Path to slope.asc raster (aspect.asc loaded if isSelectedAspect=true)
     * @return Updated costList with slope/aspect costs added
     */
    public Dictionary slopeInput(Dictionary costList, boolean isSelectedAspect, String path) throws IOException {

        Dictionary headerInfo = getHeader(path);
        double[][] slopeMatrix = getDetails(headerInfo, path);
        double[][] tempMatrixSlope = new double[slopeMatrix.length][slopeMatrix[0].length];
        double[][] aspectMatrix = new double[slopeMatrix.length][slopeMatrix[0].length];

        if (isSelectedAspect == true) {
            System.out.println("Importing Aspect Data ...");
            aspectMatrix = aspectInput();
        } else {
            for (int i = 0; i < aspectMatrix.length; i++) {
                for (int j = 0; j < aspectMatrix[0].length; j++) {
                    aspectMatrix[i][j] = -9999;
                }
            }
        }

        for (int i = 0; i < slopeMatrix.length; i++) {
            for (int j = 0; j < slopeMatrix[0].length; j++) {
                if (slopeMatrix[i][j] == (int) headerInfo.get("NoData")) {
                    tempMatrixSlope[i][j] = -9999;
                } else if (slopeMatrix[i][j] <= 0.5) {
                    tempMatrixSlope[i][j] = 0;

                } else if (slopeMatrix[i][j] > 0.5 & slopeMatrix[i][j] <= 1) {
                    tempMatrixSlope[i][j] = 0.1;

                } else if (slopeMatrix[i][j] > 1 & slopeMatrix[i][j] <= 1.5) {
                    tempMatrixSlope[i][j] = 0.2;

                } else if (slopeMatrix[i][j] > 1.5 & slopeMatrix[i][j] <= 2) {
                    tempMatrixSlope[i][j] = 0.3;

                } else if (slopeMatrix[i][j] > 2 & slopeMatrix[i][j] <= 2.5) {
                    tempMatrixSlope[i][j] = 0.4;

                } else if (slopeMatrix[i][j] > 2.5 & slopeMatrix[i][j] <= 3) {
                    tempMatrixSlope[i][j] = 0.5;

                } else if (slopeMatrix[i][j] > 3) {
                    tempMatrixSlope[i][j] = 1;
                }
            }
        }

        int indexNew = 0;
        for (int i = 0; i < tempMatrixSlope.length; i++) {
            for (int j = 0; j < tempMatrixSlope[0].length; j++) {
                indexNew = indexNew + 1;
                if (indexNew < costList.size()) {
                    double[] costs = (double[]) costList.get(indexNew);
                    double[] slopeKernel = kernel(tempMatrixSlope, i, j);
                    double[] aspectKernel = kernel(aspectMatrix, i, j);
                    costs = solveSlope(costs, slopeKernel,aspectKernel);
                    costList.put(indexNew, costs);
                }

            }
        }
        return costList;
    }


    private int[][] cellCount() throws FileNotFoundException, IOException {

        Dictionary headerInfo = getHeader("Datasets/ASCII/landcover.asc");
        int[][] cellMatrix = new int[(int)headerInfo.get("Rows")][(int)headerInfo.get("Columns")];
        int z = 1;
        for (int i = 0; i < cellMatrix.length; i++) {
            for (int j = 0; j < cellMatrix[0].length; j++) {
                cellMatrix[i][j] = z;

                z = z + 1;
            }
        }
        return cellMatrix;
    }


    public ArrayList cells() throws IOException {
        ArrayList cellList = new ArrayList();
        Dictionary headerInfo = getHeader("Datasets/ASCII/landcover.asc");
        int rows = (int)headerInfo.get("Rows");
        int columns = (int)headerInfo.get("Columns");
//        int cellAmount = rows * columns;
        int[][] cellNumber = cellCount();
//        for (int  = 0; i < 9; i++)
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                int[] cellKernel = cellKernel(cellNumber, i, j);
                cellList.add(cellKernel);
            }
        }
        return cellList;
    }

    private double adjustSlope(double slopeCost, double aspect, int zAdjust) {
        int noData = -9999;

        if (aspect == noData) return slopeCost;

        double diff = Math.abs(aspect - zAdjust);

        if (diff == 0 || diff == 4) return slopeCost + 0.2;
        if (diff == 1 || diff == 3 || diff == 5 || diff == 7) return slopeCost + 0.1;
        return slopeCost;
    }

    private double[] solveSlope(double[] costs, double[] slopeKernel, double[] aspectKernel) {

        // Center should NOT be direction-dependent if you don’t want it changed 8 times.
        // If you truly want "missing => 0", then for center just use 0 unless you have a defined center direction.
        int noData = -9999;
        if (slopeKernel[0] == noData) return costs; // Don't adjust costs to any neighbor.

        for (int j = 1; j < 9; j++) {

            if (slopeKernel[j] == noData) continue;  // Don't adjust costs[j].

            int zAdjust = (j<=4) ? (j-1)*2+1 : (j-4)*2;

            double center = adjustSlope(slopeKernel[0], aspectKernel[0], zAdjust); // Potentially add 0.2 or 0.1 to cost in slopeKernel[0]
            double neighbor = adjustSlope(slopeKernel[j], aspectKernel[j], zAdjust); // Potentially add 0.2 or 0.1 to cost in slopeKernel[j]

            double value = (center + neighbor) / 2.0; // Average costs.
            costs[j] += value;  // Add to total construction costs.
        }
        return costs;
    }

    public double [][] aspectInput() throws FileNotFoundException, IOException {

            Dictionary headerInfo = getHeader("Datasets/ASCII/aspect.asc");
            double[][] maskMatrix = getDetails(headerInfo, "Datasets/ASCII/aspect.asc");
            double[][] tempMatrix = new double[(int)headerInfo.get("Rows")][(int)headerInfo.get("Columns")];

        for (int i = 0; i < maskMatrix.length; i++) {
            for (int j = 0; j < maskMatrix[0].length; j++) {
                if (maskMatrix[i][j] == (int) headerInfo.get("NoData")) {
                    tempMatrix[i][j] = -9999;
                } else if (maskMatrix[i][j] <= 22.5) {
                    tempMatrix[i][j] = 1;

                } else if (maskMatrix[i][j] > 22.5 & maskMatrix[i][j] <= 67.5) {
                    tempMatrix[i][j] = 2;

                } else if (maskMatrix[i][j] > 67.5 & maskMatrix[i][j] <= 112.5) {
                    tempMatrix[i][j] = 3;

                } else if (maskMatrix[i][j] > 112.5 & maskMatrix[i][j] <= 157.5) {
                    tempMatrix[i][j] = 4;

                } else if (maskMatrix[i][j] > 157.5 & maskMatrix[i][j] <= 202.5) {
                    tempMatrix[i][j] = 5;

                } else if (maskMatrix[i][j] > 202.5 & maskMatrix[i][j] <= 247.5) {
                    tempMatrix[i][j] = 6;

                } else if (maskMatrix[i][j] > 247.5 & maskMatrix[i][j] <= 292.5) {
                    tempMatrix[i][j] = 7;

                } else if (maskMatrix[i][j] > 292.5 & maskMatrix[i][j] <= 337.5) {
                    tempMatrix[i][j] = 8;

                } else {
                    tempMatrix[i][j] = 1;
                }
            }
        }
        // System.out.println(Arrays.deepToString(tempMatrix));
        return tempMatrix;
    }

    public ArrayList solveConstruction(ArrayList costList, String path) throws IOException {
        Dictionary headerInfo = getHeader(path);
        double[][] landcoverMatrix = getDetails(headerInfo, "Datasets/ASCII/roads.asc");
        System.out.println("Construction Costs continue ...");
        return costList;
    }

    /**
     * Applies river crossing penalty to cost network.
     * Multiplies costs by 1.25 when route crosses a river (detected by rivers.asc raster).
     * Uses crossIncrease() to check 3x3 sub-kernels for actual feature crossing.
     *
     * @param costList Dictionary from previous aggregation step
     * @param headerInfo Header Dictionary from getHeader()
     * @param path Path to rivers.asc raster
     * @return Updated costList with river crossing penalties applied
     */
    public Dictionary addRiverCrossings(Dictionary costList, Dictionary headerInfo, String path) throws IOException {
        Dictionary roadInfo = getHeader(path);
        double[][] roadMatrix = getDetails(roadInfo, path);
        double[][] tempMatrix = new double[(int)headerInfo.get("Rows")][(int)headerInfo.get("Columns")];
        int indexNew = 0;
        double weight = 1.25;
        for (int i = 0; i < tempMatrix.length; i++) {
            for (int j = 0; j < tempMatrix[0].length; j++) {
                indexNew = indexNew + 1;
                if (indexNew < costList.size()) {
                    double[] costs = (double[]) costList.get(indexNew);
                    double[] newCosts = crossIncrease(costs,roadMatrix, i, j, weight);
                    costList.put(indexNew, newCosts);
                }

            }
        }
        return costList;
    }

    /**
     * Applies road crossing penalty to cost network.
     * Multiplies costs by 1.25 when route crosses a road (detected by roads.asc raster).
     * Uses crossIncrease() to check 3x3 sub-kernels for actual feature crossing.
     *
     * @param costList Dictionary from previous aggregation step
     * @param headerInfo Header Dictionary from getHeader()
     * @param path Path to roads.asc raster
     * @return Updated costList with road crossing penalties applied
     */
    public Dictionary addRoadCrossings(Dictionary costList, Dictionary headerInfo, String path) throws IOException {
        Dictionary roadInfo = getHeader(path);
        double[][] roadMatrix = getDetails(roadInfo, path);
        double[][] tempMatrix = new double[(int)headerInfo.get("Rows")][(int)headerInfo.get("Columns")];
        int indexNew = 0;
        double weight = 1.25;
        for (int i = 0; i < tempMatrix.length; i++) {
            for (int j = 0; j < tempMatrix[0].length; j++) {
                indexNew = indexNew + 1;
                if (indexNew < costList.size()) {
                    double[] costs = (double[]) costList.get(indexNew);
                    double[] newCosts = crossIncrease(costs,roadMatrix, i, j, weight);
                    costList.put(indexNew, newCosts);
                }

            }
        }
        return costList;
    }

    /**
     * Applies railroad crossing penalty to cost network.
     * Multiplies costs by 1.25 when route crosses a railroad (detected by railroads.asc raster).
     * Uses crossIncrease() to check 3x3 sub-kernels for actual feature crossing.
     *
     * @param costList Dictionary from previous aggregation step
     * @param headerInfo Header Dictionary from getHeader()
     * @param path Path to railroads.asc raster
     * @return Updated costList with railroad crossing penalties applied
     */
    public Dictionary addRailCrossings(Dictionary costList, Dictionary headerInfo, String path) throws IOException {
        Dictionary roadInfo = getHeader(path);
        double[][] roadMatrix = getDetails(roadInfo, path);
        double[][] tempMatrix = new double[(int)headerInfo.get("Rows")][(int)headerInfo.get("Columns")];
        int indexNew = 0;
        double weight = 1.25;
        for (int i = 0; i < tempMatrix.length; i++) {
            for (int j = 0; j < tempMatrix[0].length; j++) {
                indexNew = indexNew + 1;
                if (indexNew < costList.size()) {
                    double[] costs = (double[]) costList.get(indexNew);
                    double[] newCosts = crossIncrease(costs,roadMatrix, i, j, weight);
                    costList.put(indexNew, newCosts);
                }

            }
        }
        return costList;
    }

    /**
     * Applies pipeline corridor DISCOUNT to cost network.
     * Unlike infrastructure crossings (penalties), existing pipelines offer cost reduction (0.75x multiplier).
     * Uses rowDecrease() with 3x3 sub-kernel alignment detection to apply discount only along
     * pipeline-aligned paths. This incentivizes route placement within existing right-of-way.
     *
     * @param costList Dictionary from previous aggregation step
     * @param headerInfo Header Dictionary from getHeader()
     * @param path Path to pipelines.asc raster
     * @return Updated costList with pipeline corridor discounts applied
     */
    public Dictionary addPipelineCorridor(Dictionary costList, Dictionary headerInfo, String path) throws IOException {
        Dictionary pipelineInfo = getHeader(path);
        double[][] pipeMatrix = getDetails(pipelineInfo, path);
        int indexNew = 0;
        double weight = .75;
        for (int i = 0; i < pipeMatrix.length; i++) {
            for (int j = 0; j < pipeMatrix[0].length; j++) {
                indexNew = indexNew + 1;
                if (indexNew < costList.size()) {
                    double[] costs = (double[]) costList.get(indexNew);
                    double[] newCosts = rowDecrease(pipeMatrix,costs, i, j, weight);
                    costList.put(indexNew, newCosts);
                }

            }
        }
        return costList;
    }



    public double[] rowDecrease(double[][] matrix, double[]costs, int i, int j, double weight) {


        int z = (i * 3) + 1;
        int d = (j * 3) + 1;

        double[] kernel = kernel(matrix, z, d);

        for (int r = 1; r < 9; r++) {


            double[] kernel4 = kernel(matrix, (z + rows[r] * 3), (d + cols[r] * 3));

            //case 1 slide 2
            if ((kernel[1] > 0 && kernel4[3] > 0) || (kernel[3] < 0 && kernel4[1] < 0)) {
                costs[r] = costs[r] * weight;

            } else if ((kernel[2] > 0 && kernel4[4] > 0) || (kernel[4] < 0 && kernel4[2] < 0)) {
                costs[r] = costs[r] * weight;

            } else if ((kernel[6] > 0 && kernel4[8] > 0) || (kernel[8] < 0 && kernel4[6] < 0)) {
                costs[r] = costs[r] * weight;

            } else if ((kernel[5] > 0 && kernel4[7] > 0) || (kernel[7] < 0 && kernel4[5] < 0)) {
                costs[r] = costs[r] * weight;

            } else if ((kernel[5] > 0 && kernel4[6] > 0) || (kernel[6] < 0 && kernel4[5] < 0)) {
                costs[r] = costs[r] * weight;

            } else if ((kernel[7] > 0 && kernel4[8] > 0) || (kernel[8] < 0 && kernel4[7] < 0)) {
                costs[r] = costs[r] * weight;

            } else {
                costs[r] = costs[r] * 1;
            }

        }
        return  costs;
    }

    public double[] crossIncrease(double[] costs, double[][] matrix, int i, int j, double weight) {

        double[] increaseValue = new double[9];
        int z = (i * 3) + 1;
        int d = (j * 3) + 1;
        double[] kernel = kernel(matrix, z, d);

        for (int r = 1; r < 9; r++) {

            switch(r){

            case 1:
                double[] kernel2 = kernel(matrix, z - 3, d);
                //Case 1
                if ((kernel[0] < 0 && kernel[1] < 0 && kernel2[0] < 0 && kernel2[3] < 0)) {

                }
                //Case 2
                else if ((kernel[0] < 0 && kernel[2] < 0 && kernel[6] < 0 && kernel2[0] < 0 && kernel2[2] < 0 && kernel[7] < 0)) {

                }
                //Case 3
                else if ((kernel[0] < 0 && kernel[4] < 0 && kernel[5] < 0 && kernel2[0] < 0 && kernel2[4] < 0 && kernel[8] < 0)) {

                }
                //Case 4
                else if ((kernel[0] < 0 && kernel[2] < 0 && kernel[6] < 0 && kernel2[0] < 0 && kernel2[3] < 0 && kernel[7] < 0)) {

                }
                //Case 5
                else if ((kernel[0] < 0 && kernel[4] < 0 && kernel[5] < 0 && kernel2[0] < 0 && kernel2[3] < 0 && kernel[8] < 0)) {

                } else {
                    costs[r] = costs[r] * weight;
                }
            break;

            case 2:
                double[] kernel3 = kernel(matrix, z, d + 3);
                //Case 1
                if ((kernel[0] < 0 && kernel[2] < 0 && kernel[0] < 0 && kernel3[4] < 0)) {

                }
                //Case 2
                else if ((kernel[0] < 0 && kernel[3] < 0 && kernel[7] < 0 && kernel3[0] < 0 && kernel3[3] < 0 && kernel[8] < 0)) {

                }
                //Case 3
                else if ((kernel[0] < 0 && kernel[1] < 0 && kernel[6] < 0 && kernel3[0] < 0 && kernel3[1] < 0 && kernel[5] < 0)) {

                }

                //Case4
                else if ((kernel[0] < 0 && kernel[1] < 0 && kernel[6] < 0 && kernel3[0] < 0 && kernel3[4] < 0 && kernel[5] < 0)) {

                }
                //Case5
                else if ((kernel[0] < 0 && kernel[3] < 0 && kernel[7] < 0 && kernel3[0] < 0 && kernel3[4] < 0 && kernel[8] < 0)) {

                } else {
                    costs[r] = costs[r] * weight;
                }

            break;

            case 3:
                double[] kernel4 = kernel(matrix, z + 3, d);
                //Case 1
                if ((kernel[0] < 0 && kernel[3] < 0 && kernel4[0] < 0 && kernel4[1] < 0)) {

                }
                //Case 2
                else if ((kernel[0] < 0 && kernel[2] < 0 && kernel[7] < 0 && kernel4[0] < 0 && kernel4[2] < 0 && kernel[6] < 0)) {

                }
                //Case 3
                else if ((kernel[0] < 0 && kernel[4] < 0 && kernel[8] < 0 && kernel4[0] < 0 && kernel4[4] < 0 && kernel[5] < 0)) {

                }
                //Case 4
                else if ((kernel[0] < 0 && kernel[2] < 0 && kernel[6] < 0 && kernel4[0] < 0 && kernel4[1] < 0 && kernel[7] < 0)) {

                }
                //Case 5
                else if ((kernel[0] < 0 && kernel[4] < 0 && kernel[8] < 0 && kernel4[0] < 0 && kernel4[1] < 0 && kernel[5] < 0)) {

                } else {
                    costs[r] = costs[r] * weight;
                }
            break;

            case 4:
                  double[] kernel5 = kernel(matrix, z, d - 3);
                  //Case 1
                  if ((kernel[0] < 0 && kernel[2] < 0 && kernel5[0] < 0 && kernel5[4] < 0)) {

                  }
                  //Case 2
                  else if ((kernel[0] < 0 && kernel[3] < 0 && kernel[8] < 0 && kernel5[0] < 0 && kernel5[3] < 0 && kernel[7] < 0)) {

                  }
                  //Case 3
                  else if ((kernel[0] < 0 && kernel[1] < 0 && kernel[5] < 0 && kernel5[0] < 0 && kernel5[1] < 0 && kernel[6] < 0)) {

                  }

                  //Case4
                  else if ((kernel[0] < 0 && kernel[1] < 0 && kernel[5] < 0 && kernel5[0] < 0 && kernel5[2] < 0 && kernel[6] < 0)) {

                  }
                  //Case5
                  else if ((kernel[0] < 0 && kernel[3] < 0 && kernel[3] < 0 && kernel5[0] < 0 && kernel5[7] < 0 && kernel[3] < 0)) {

                  } else {
                      costs[r] = costs[r] * weight;
                  }
            break;
//
            case 5:
                  double[] kernel6 = kernel(matrix, z, d - 3);
                  //Case 1
                  if ((kernel[0] < 0 && kernel[5] < 0 && kernel6[0] < 0 && kernel6[7] < 0)) {

                  } else {
                      costs[r] = costs[r] * weight;
                  }

            break;

            case 6:
                double[] kernel7 = kernel(matrix, z, d - 3);
                //Case 1
                if ((kernel[0] < 0 && kernel[6] < 0 && kernel7[8] < 0 && kernel7[0] < 0)) {

                } else {
                    costs[r] = costs[r] * weight;
                }

            break;
            case 7:
                double[] kernel8 = kernel(matrix, z, d - 3);
                //Case 1
                if ((kernel[0] < 0 && kernel[7] < 0 && kernel8[5] < 0 && kernel8[0] < 0)) {

                } else {
                    costs[r] = costs[r] * weight;
                }

            break;
            case 8:
                double[] kernel9 = kernel(matrix, z, d - 3);
                //Case 1
                if ((kernel[0] < 0 && kernel[8] < 0 && kernel9[6] < 0 && kernel9[0] < 0)) {

                } else {
                    costs[r] = costs[r] * weight;
                }
            break;
            }
        }

        return costs;
    }

    /**
     * Computes landcover costs between center cell and each neighbor.
     * Averages the center cell cost with each neighbor's cost: (center + neighbor) / 2.
     * Preserves -9999 (no-data) values for boundary cells and no-data cells.
     *
     * @param costKernel 9-element array from kernel() with landcover base weights
     * @return double[9] array with averaged costs; costs[0]=center averaged with itself (== center value)
     */
    public double[] solveLand(double[] costKernel) {
        final double noData = -9999;
        int n = costKernel.length;
        double[] costs = new double[n];

        for (int j = 0; j < n; j++) {
            if (costKernel[j] == noData || costKernel[0] == noData) {
                costs[j] = noData;
            } else {
                costs[j] = 0.5 * (costKernel[0] + costKernel[j]);
            }
        }

        return costs;
    }




    /**
     * Extracts a 9-cell kernel neighborhood around a grid cell.
     * Returns the center cell and its 8 neighbors in a fixed array layout.
     *
     * Kernel layout (indices 0-8):
     * <pre>
     *   [5] [1] [6]
     *   [4] [0] [2]
     *   [8] [3] [7]
     * </pre>
     * Where: 0=center, 1=north, 2=east, 3=south, 4=west, 5=NW, 6=NE, 7=SE, 8=SW
     *
     * @param array 2D raster grid
     * @param i Row index (0-based, from top)
     * @param j Column index (0-based, from left)
     * @return double[9] array with cell values; boundary cells return -9999 for out-of-bounds neighbors
     */
    public double[] kernel(double[][] array, int i, int j) {

        double noData = -9999;
        double[] kernel = new double[9];

        for (int r = 0; r < 9; r++) {
            int ii = i + rows[r];
            int jj = j + cols[r];

            if (ii < 0 || ii >= array.length || jj < 0 || jj >= array[0].length) {
                kernel[r] = noData;
            } else {
                kernel[r] = array[ii][jj];
            }
        }

        return kernel;
    }

    // Kernel for calculations
    public int[] cellKernel(int[][] array, int i, int j) {

        int noData = -9999;
        int[] kernelInt = new int[9];

        for (int r = 0; r < 9; r++) {
            int row = rows[r];
            int col = cols[r];

            int ii = i + row;
            int jj = j + col;

            if (ii < 0 || ii >= array.length || jj < 0 || jj >= array[0].length) {
                kernelInt[r] = noData;
            } else {
                kernelInt[r] = array[ii][jj];
            }
        }

        return kernelInt;
    }

    private int getActiveNodes(Dictionary headerInfo, double[][] matrix) {
        int activeNodes = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (matrix[i][j] != (int) headerInfo.get("NoData")) {
                    activeNodes++;
                }
            }
        }

        return activeNodes;

    }

    public double round(double value, int places) {
        if (places < 0) {
            throw new IllegalArgumentException();
        }
        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }


    /**
     * Exports cost network to tab-delimited text file for SimCCS pipeline routing.
     * Writes grid metadata header followed by node-cost pairs.
     * Filters out negative/no-data costs (non-traversable paths).
     * Only includes neighbor connections with positive cost (>0).
     *
     * Output format:
     * - Header: All Nodes, Active Nodes, nCols, nRows, xllCorner, yllCorner, cellSize, NODATA_value
     * - Nodes: Tab-separated cell IDs for current cell and reachable neighbors
     * - Costs: Tab-separated costs corresponding to neighbor cells
     * - Example line: 37737\\t37517\\t37738\\t37957 (cells) -> \\t2.15\\t3.42\\t2.08\\t3.45 (costs)
     *
     * @param costList Dictionary from final aggregation pipeline
     * @param headerInfo Header Dictionary from getHeader()
     * @param outPut BufferedWriter to output file (typically Outputs/Construction Costs.txt)
     * @throws IOException if file write fails
     */
    public void writeTxt(Dictionary costList, Dictionary headerInfo, BufferedWriter outPut) throws IOException {

        double[][] cellMatrix = getDetails(headerInfo, "Datasets/ASCII/landcover.asc");
        outPut.write("All Nodes" + "\t" + ((int) headerInfo.get("Rows") * (int) headerInfo.get("Columns")));
        outPut.newLine();
        outPut.write("Active Nodes" + "\t" + getActiveNodes(headerInfo, cellMatrix));
        outPut.newLine();
        outPut.write("nCols" + "\t" + headerInfo.get("Columns"));
        outPut.newLine();
        outPut.write("nRows" + "\t" + headerInfo.get("Rows"));
        outPut.newLine();
        outPut.write("xllCorner" + "\t" + headerInfo.get("xllCorner"));
        outPut.newLine();
        outPut.write("yllCorner" + "\t" + headerInfo.get("yllCorner"));
        outPut.newLine();
        outPut.write("cellSize" + "\t" + headerInfo.get("CellSize"));
        outPut.newLine();
        outPut.write("NODATA_value" + "\t" + headerInfo.get("NoData"));
        outPut.newLine();
        int rows = (int) headerInfo.get("Rows");
        int cols = (int) headerInfo.get("Columns");

        for (int i = 1; i < costList.size(); i++) {

            double[] costs = (double[]) costList.get(i);
            int[] cells = new int[9];
            if (costs[0] != (int) headerInfo.get("NoData")) {

                cells[0] = i;
                cells[1] = i - cols;
                cells[2] = i + 1;
                cells[3] = i + cols;
                cells[4] = i - 1;
                cells[5] = i - (cols + 1);
                cells[6] = i - (cols - 1);
                cells[7] = i + (cols + 1);
                cells[8] = i + (cols - 1);

//                System.out.println(Arrays.toString(cells));
//                System.out.println(Arrays.toString(costs));

                //Prints the cells
                ArrayList<Integer> printCells = new ArrayList<Integer>();
                printCells.add(cells[0]);
                for (int z = 1; z < 9; z++) {
                    if (costs[z] > 0) {
                        printCells.add(cells[z]);

                    }}
                    // printing the cells
                    if (printCells.size() > 1) {
                        String x = printCells.toString().replace('[', ' ');
                        x = x.replace(']', ' ');
                        x = x.trim().replace(',', '\t');
                        outPut.write(x);
                        outPut.newLine();
                    }
                }

//
//            //Prints the costs
                List<Double> printCosts = new ArrayList<>();

                for (int z = 1; z < 9; z++) {

                    if (costs[z] > 0) {
                        printCosts.add(round(costs[z], 2));

                    }
                }
//
                if (printCosts.size() > 0) {
                    String x = printCosts.toString().replace('[', ' ');

                    x = x.replace(']', ' ');
                    x = x.replace(',', '\t');
                    outPut.write("\t");
                    outPut.write(x);
                    outPut.newLine();
                }
            }


        outPut.close();
    }

    /**
     * Computes per-row distance multipliers for each compass direction.
     * Calculates actual geographic cell dimensions using haversine formula,
     * accounting for latitude effects on longitude distance.
     * Returns a matrix used by solveDistance() to scale costs geographically.
     *
     * @param headerInfo Header Dictionary from getHeader() containing coordinates and cellSize
     * @return double[rows][4] where each row contains multipliers for:
     *         [i][0]=E/W distance, [i][1]=N/S distance, [i][2]=diagonal NE/SW, [i][3]=diagonal NW/SE
     */
    public double[][] distanceMultiplier(Dictionary headerInfo) {

        int rows = (int) headerInfo.get("Rows");
        double cellSize = (double) headerInfo.get("CellSize");
        double yllCorner = (double) headerInfo.get("yllCorner");
        double xllCorner = (double) headerInfo.get("xllCorner");

        double[][] cellMatrix = new double[rows][4];

        //x multiplier
        for (int i = 0; i < rows; i++) {

            double lat1 = yllCorner + ((rows - (i + 1)) * cellSize) + (cellSize / 2);
            double lat2 = yllCorner + ((rows - (i + 1)) * cellSize) + (cellSize / 2);
            double lon1 = xllCorner;
            double lon2 = xllCorner + cellSize;
            cellMatrix[i][0] = haversineDistance(lat1, lon1, lat2, lon2);

        }

        //y multiplier
        for (int i = 0; i < rows; i++) {

            double lat1 = yllCorner;
            double lat2 = yllCorner + cellSize;
            double lon1 = xllCorner;
            double lon2 = xllCorner;

            cellMatrix[i][1] = haversineDistance(lat1, lon1, lat2, lon2);
        }

        //xy-up
        for (int i = 0; i < rows; i++) {

            double lat1 = yllCorner + ((rows - (i + 1)) * cellSize) + (cellSize / 2);
            double lat2 = (yllCorner + ((rows - (i + 1)) * cellSize) + (cellSize / 2)) + cellSize;
            double lon1 = xllCorner;
            double lon2 = xllCorner + cellSize;

            cellMatrix[i][2] = haversineDistance(lat1, lon1, lat2, lon2);

        }

        //XY-down
        for (int i = 0; i < rows; i++) {

            double lat1 = yllCorner + ((rows - (i + 1)) * cellSize) + (cellSize / 2);
            double lat2 = (yllCorner + ((rows - (i + 1)) * cellSize) + (cellSize / 2)) - cellSize;
            double lon1 = xllCorner;
            double lon2 = xllCorner + cellSize;

            cellMatrix[i][3] = haversineDistance(lat1, lon1, lat2, lon2);

        }
        return cellMatrix;
    }

    /**
     * Aggregates the 9-cell cost kernel into a single per-pixel value.
     * Computes the average cost to traverse to valid (non-negative) neighbors,
     * excluding the center cell and any invalid cells marked as -9999.
     *
     * @param cellKernel double[9] array where index 0=center, indices 1-8=neighbors
     * @return Average cost to valid neighbors; 0 if no valid neighbors exist
     */
    public double aggregateCellCost(double[] cellKernel) {
        double sum = 0;
        int validCount = 0;

        // Skip index 0 (center), iterate through neighbors (1-8)
        for (int i = 1; i < cellKernel.length; i++) {
            if (cellKernel[i] >= 0) {  // Exclude no-data/invalid values
                sum += cellKernel[i];
                validCount++;
            }
        }

        return validCount > 0 ? sum / validCount : 0;
    }

}
