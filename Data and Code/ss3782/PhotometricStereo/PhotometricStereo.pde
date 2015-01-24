///
/// CS-583: Introduction to Computer Vision
/// https://www.cs.drexel.edu/~kon/introcompvis/
///
/// Project 3 - Photometric Stereo
///
/// @note Look for the tag 'TODO' to see what you must complete
///

// Needed for reading files
import java.io.FileNotFoundException;
import java.util.*;
// Needed for the various Matrix and Vector operations
import no.uib.cipr.matrix.*;
import no.uib.cipr.matrix.sparse.*;
// Needed for functions like Max and Min
import java.lang.Math.*;
// Needed for the .ini file parsing
import org.ini4j.*;
import java.util.Formatter;




//
// Constants
//

// ini file reading constants including filename and key names
// static final String DEFAULT_INI_FILE = "panorama.ini";

// Homographic indexes
static final int X = 0;
static final int Y = 1;
static final int Z = 2;

// Whether or not to print debug information
static final boolean DEBUG = true;

// Size of red dot to draw. Used in showPoints
static final int DOT_SIZE = 10;

// ini file reading constants including filename and key names
static final String DEFAULT_INI_FILE = "photometricStereo.ini";
static final String CHROME = "Chrome";
static final String PATH = "path";
static final String OUTPUT = "output";
static final String SUBJECTS = "Subjects";
static final String MASK = ".mask.png";
static final String PATTERN = ".%d.png";

// System constants
static final String SLASH = System.getProperty("file.separator");
static final int MAX_IMAGES = 12;

//
static final boolean DO_CHROME =  true;
static final boolean DO_ALBEDOS = true;
static final boolean DO_NORMALS = true;
static final boolean DO_DEPTH =   true;

// Globals for the images to show
PImage imageToShow = createImage(400, 400, RGB);;
Matrix pointsToShow = new DenseMatrix(0,0);

//
// Functions
//


///
/// Set up screen size, turn off re-drawing loop
///
void setup() {
  size(400, 400);
  noLoop();

  // Read .ini file
  Ini ini = new Ini();
  try {
    BufferedReader reader = createReader(DEFAULT_INI_FILE);
    if(reader == null) {
      throw new FileNotFoundException(DEFAULT_INI_FILE);
    }
    ini.load(reader);

    if(DO_CHROME) {
      //
      // TODO : Estimate light direction for each chrome image using 'computeOneLightDirection'
      // Write each direction to a file (suggestion, use filename in ini file)
      // (Hint: use Processing's 'createWriter' to make a PrintWriter)
      //
      println("Analyzing chrome images to determine light sources");
      Ini.Section section = ini.get(CHROME);
      PImage[] chromeImages = loadImageDirectory(section.get(PATH), section.get(PATH) + PATTERN);
      PImage chromeMask = loadImage(section.get(PATH) + SLASH + section.get(PATH) + MASK);
      
      PrintWriter Writelights = createWriter(section.get(OUTPUT));
      //for (PImage imgs : chromeImages)
     for(int i=0;i<chromeImages.length;i++) 
      {
        
        DenseVector writing = computeOneLightDirection(chromeImages[i], chromeMask);
        Writelights.println(writing.get(X) + " " + writing.get(Y) + " " + writing.get(Z));
      }
      Writelights.close();

      
//DenseVector light = new DenseVector(2);
      //for(int i=0; i<chromeImages.length;i++)
      
//      light = computeOneLightDirection(chromeImages, chromeMask);
      

    } else {
      println("Skipping light direction calculation section");
    }

    //
    // Iterate over each subject and perform the requested actions
    //
    Iterator sectionIterator = ini.get(SUBJECTS).getAll(PATH).iterator();
    while(sectionIterator.hasNext()) {
      String path = (String)sectionIterator.next();
      println("Now working with " + path + " images");
      
      DenseMatrix lights=readPoints3D("lights.txt");
      PImage[] images=loadImageDirectory(path, path+PATTERN);
      PImage mask=loadImage(path + SLASH + path + MASK);
      
      if(DO_NORMALS) {
        //
        // TODO: Compute normal map image by calling 'computeNormals'. Save the result.
        //
        println("Computing normal map");
    PImage normals = computeNormals(images, mask, lights);
    normals.save(path + "_normal.png");

      } else {
        println("Skipping normal calculation section");
      }
      if(DO_ALBEDOS) {
        //
        // TODO: Compute albedo image by calling 'computeAlbedos'. Save the result.
        //
        println("Computing albedo map");
        PImage albedos = computeAlbedos(images, mask, lights);
        albedos.save(path + "_albedo.png");



      } else {
        println("Skipping albedo calculation section");
      }
      if(DO_DEPTH) {
        //
        // TODO: Compute depth image by calling 'computeDepths'. Save the result.
        //
        println("Computing depth map");
        PImage depths = computeDepths(images, mask, lights);
        depths.save(path + "_depth.png");

      }

      println("Completed " + path + " subject.");
    }

  } catch(FileNotFoundException e) {
    println("You must have the file '" + e.getMessage() + "' present.");
  } catch(IOException e) {
    println("Failure reading ini file: " + e.getMessage());
    exit();
  }
  println("Done.");
}

///
/// Draw function
///
void draw() {
  background(0);

  float divideBy = 1.0;
  Matrix newPoints = pointsToShow.copy();
  if(imageToShow.width > width || imageToShow.height > height) {
    divideBy = max(float(imageToShow.width)/width, float(imageToShow.height)/height);
    newPoints.scale(1.0/divideBy);
  }

  int originX = floor((width-imageToShow.width/divideBy)/2);
  int originY = floor((height-imageToShow.height/divideBy)/2);

  for(int r = 0; r < pointsToShow.numRows(); r++) {
    newPoints.set(r, X, newPoints.get(r, X) + originX);
    newPoints.set(r, Y, newPoints.get(r, Y) + originY);
  }

  // Show the picture scaled appropriately
  image(imageToShow, originX, originY, imageToShow.width/divideBy, imageToShow.height/divideBy);

  // Show the points
  showPoints(newPoints);
}


///
/// Given an image of a chrome sphere under a point light source, and an
/// associated mask, returns the direction of the point light source.
///
/// @note assumes the viewer is at (0, 0, 1)
/// @param image image of a chrome sphere under a point light source
/// @param mask mask that indicates where the sphere is in the image
///
DenseVector computeOneLightDirection(PImage image, PImage mask) {
  // Find the center and radius of the sphere from the mask
  DenseVector sphereCenter = new DenseVector(2);
  float sphereRadius = findCenterAndRadius(mask, sphereCenter);
  println("RADIUS" + sphereRadius);

  // TODO: Compute a weighted average of the brightest pixel locations.
  // (You might want to set the threshold for finding the highlight region to 0.9.)
  DenseVector highlightCenter = new DenseVector(new double[]{0.0, 0.0});

  double xs=0;
  for(int i=0;i<image.height;i++)
  {
    for(int j=0;j<image.width;j++)
    {
//double x=brightness(image.get(i,j));

//double x=red(image.get(i,j));
  double x=bright(image.get(i,j));
      if(x>0.9)
      {
        highlightCenter.add(X, i*x);
        highlightCenter.add(Y, j*x);
        xs+=x;
      }
    }
  }
  highlightCenter.scale(1/xs);
  
  println("HIGHLIGHTCENTER:"+ highlightCenter);

  // TODO: Compute the surface normal at the point
  DenseVector normal = normalOnSphere(highlightCenter, sphereCenter, sphereRadius);

  // TODO: Compute light source direction from normal
  DenseVector lightSource = new DenseVector(3);
  // s = 2(n.R)n - R
  //reflection Direction R=[0 0 1.0]
  DenseVector R = new DenseVector(new double[]{0.0, 0.0, 1.0});
  lightSource.set(normal);
  lightSource.scale(2*(normal.dot(R)));
 // lightSource.add(normal, -1.0*R);
  lightSource.add(Z, -1.0);
  
 
println("LIGHTSOURCE" + lightSource);

  return lightSource;
}


///
/// Assuming the 'image' is a mask that has the shape of a sphere, sets the center parameters to
/// and returns the radius (in pixels)
///
/// @param mask mask of a sphere
/// @param center set by this method to the center of the sphere
/// @post center is modified to contain the proper information
/// @returns the radius of the sphere
///
float findCenterAndRadius(PImage mask, DenseVector center) {
  //
  // TODO: Compute 'center' and return radius
  //
  float radius = 0;
// center={0.0,0.0};
  center.set(X,0);
  center.set(Y,0);
  double all=0;
  
  int minRow=mask.width, minCol=mask.height, maxRow=0, maxCol=0;
 
 
  // minRow, maxRow, minCOl, maxCol calculations
  for(int j=0;j<mask.height;j++)
  {
    for(int i=0;i<mask.width;i++)
    {
      int z=mask.get(i,j);
      double s=red(z);
  center.add(X, s*i);
  center.add(Y, s*j);
    all+=s;
   if(s>0)
   {
     minRow=min(i,minRow);
     minCol=min(j,minCol);
     maxRow=max(i,maxRow);
     maxCol=max(j,maxCol);
   }
    }
  }
  center.scale(1/all);
  println("minRow:" + minRow);
  println("minCol:" + minCol);
  println("maxRow:" + maxRow);
  println("maxCol:" + maxCol);
  
  center.set(X, (minRow+maxRow)/2);
  center.set(Y, (minCol+maxCol)/2);
  println("Center:" + center);
  
  radius=(maxRow-minRow)/2;
 // println(" Radius"+radius);
  return radius;
  
}



double bright(int a)
{
  double RED,GREEN,BLUE;
  RED=red(a)/255.0;
  GREEN=green(a)/255.0;
  BLUE=blue(a)/255.0;
  double result = sqrt((float)RED*(float)RED+(float)GREEN*(float)GREEN+(float)BLUE*(float)BLUE);
  return result;
}
///
/// Given the size and location of a sphere, and the location of a highlight
/// computes the normal of the sphere at that point
///
/// @param highlightCenter coordinates of the center of the highlight
/// @param sphereCenter coordinates of the center of the sphere
/// @param radius radius of the sphere
/// @returns a 3D vector indicating the surface normal of the sphere at the point of the highlight
///
DenseVector normalOnSphere(DenseVector highlightCenter, DenseVector sphereCenter, float radius) {
  //
  // TODO: Compute normal
  //
  DenseVector normal = new DenseVector(3); 
  //double rad=radius*1.0;
  float px=(float)highlightCenter.get(X) - (float)sphereCenter.get(X);
  float py=((float)highlightCenter.get(Y) - (float)sphereCenter.get(Y));
  double z= sqrt(sq(radius)-sq(px)-sq(py));
  normal.set(X, px); 
  normal.set(Y, py);
  normal.set(Z, z);
  normal.scale(1/radius);
  println("NORMAL" + normal);
  return normal;
}


///
/// Computes the normal at each pixel and constructs an image visualizing them
///
/// @param images array of images lit from different directions
/// @param mask mask indicating where the subject is
/// @param lights the light directions corresponding to each image
/// @returns an image where each pixel of the subject has been colored based on the surface normal
/// @pre 'images' has as many entries as 'lights' has rows
/// @pre 'lights' has 3 columns
/// @pre 'mask' and each entry in 'images' are all the same size
///
PImage computeNormals(PImage[] images, PImage mask, DenseMatrix lights) {
  PImage result = createImage(mask.width, mask.height, RGB);
  for(int y = 0; y < mask.height; y++) {
    for(int x = 0; x < mask.width; x++) {
      if(red(mask.get(x, y)) > 0) {
        DenseVector normal = computeOneNormal(images, lights, x, y);
        normal.add(new DenseVector(new double[]{1,1,1}));
        normal.scale(255 * 0.5);
        result.set(x, y, color(Math.round(normal.get(X)), Math.round(normal.get(Y)), Math.round(normal.get(Z))));
      }
    }
  }
  return result;
}


///
/// Given the pixel values and light directions computes the normal to the surface
/// (assumes Lambertian reflection)
///
/// @param images array of images of a subject lit from different directions
/// @param lights light sources corresponding to each image
/// @param x horizontal component pixel where normal should be calculated
/// @param y vertical  component pixel where normal should be calculated
/// @returns a 3-dimensional vector which is the surface normal at (x,y)
///
DenseVector computeOneNormal(PImage[] images,  DenseMatrix lights, int x, int y) {
  //
  // TODO: Compute the surface normal at pixel (x,y)
  // (Try the weighted least-squares explained in the class for handling shadows.)
  //
  DenseVector normal = new DenseVector(3);
  int len=images.length;
  
  // We have to find the pseudoinverse: ~n=(pow((StS),-1))*St*I
  DenseVector I=new DenseVector(len);
  DenseMatrix S= lights.copy();
  for(int i=0;i<len;i++)
  {
    //    double intensity=brightness(images[i].get(x,y));

    double intensity=bright(images[i].get(x,y));
    I.set(i, intensity*intensity);
    
    for(int j=0;j<3;j++)
    {
      double var=S.get(i,j);
      S.set(i,j, intensity*var);
    }
  }
  
  DenseMatrix St=new DenseMatrix(S.numColumns(), S.numRows());
  
  S.transpose(St);
  
  DenseMatrix StSinv=new DenseMatrix(S.numColumns(), S.numRows());
  
  DenseMatrix StS=new DenseMatrix(S.numColumns(), S.numColumns());
  St.mult(S, StS);
  
  try {
  StS.solve(St, StSinv);
  StSinv.mult(I, normal);
  //DenseVector normalisation=sqrt(add(X,normal*normal));
    //normal.scale(1/normalisation);

  normal.scale(1/normal.norm(no.uib.cipr.matrix.Vector.Norm.TwoRobust)); //sqrt(sum(sqrs))
  //normal.scale(1/normalisation);
  } catch(MatrixSingularException exception) {
    normal.zero();
  }
 
     //println("NORMALLY :" + normal);
   return normal;

}


///
/// Computes the albedos at each pixel and constructs an image visualizing them
///
/// @param images array of images lit from different directions
/// @param mask mask indicating where the subject is
/// @param lights the light directions corresponding to each image
/// @returns an image where each pixel of the subject has been colored based on the surface albedo
/// @pre 'images' has as many entries as 'lights' has rows
/// @pre 'lights' has 3 columns
/// @pre 'mask' and each entry in 'images' are all the same size
///
PImage computeAlbedos(PImage[] images, PImage mask, DenseMatrix lights) {
  PImage result = createImage(mask.width, mask.height, RGB);
  for(int y = 0; y < mask.height; y++) {
    for(int x = 0; x < mask.width; x++) {
      if(red(mask.get(x, y)) > 0) {
        DenseVector normal = computeOneNormal(images, lights, x, y);
        DenseVector albedo = computeOneAlbedo(images, lights, x, y, normal);
        // Hack to scale
        albedo.scale(0.75 * 255);
        result.set(x, y, color(Math.round(albedo.get(X)), Math.round(albedo.get(Y)), Math.round(albedo.get(Z))));
      }
    }
  }
  return result;
}


///
/// Given the pixel values, light directions, and normal computes the albedo for a given location
/// (assumes Lambertian reflection)
///
/// @param images array of images of a subject lit from different directions
/// @param lights light sources corresponding to each image
/// @param x horizontal component pixel where normal should be calculated
/// @param y vertical component pixel where normal should be calculated
/// @param normal surface nurmal at x,y
/// @returns a 3-dimensional vector which is the surface normal at (x,y)
///
DenseVector computeOneAlbedo(PImage[] images, DenseMatrix lights, int x, int y, DenseVector normal) {
  //
  // TODO: Compute the surface normal at pixel (x,y)
  // (Try the weighted least-squares explained in the class for handling shadows.)
  // We need to find: rho=transpose[rhoR, rhoG, rhoB]
  DenseVector rho = new DenseVector(3);
  int len=images.length;
  DenseVector SurfNorm = new DenseVector(len);
  lights.mult(normal, SurfNorm);
  
  // Xr = I
 DenseMatrix X= new DenseMatrix(len*3, 3);
 DenseVector I= new DenseVector(len*3);
 
 for(int j=0;j<len;j++)
 {
   double Surf_j= SurfNorm.get(j);
   int a =3*j;
   int clr = images[j].get(x,y);
    X.set(a+0, 0, Surf_j);
    X.set(a+1, 1, Surf_j);
    X.set(a+2, 2, Surf_j);
    I.set(a+0, red(clr)/255.0);
    I.set(a+1, green(clr)/255.0);
    I.set(a+2, blue(clr)/255.0);
 }
 
 X.solve(I,rho);
 
  return rho;
}


///
/// Computes the depths at each pixel and constructs an image visualizing them
///
/// @param images array of images lit from different directions
/// @param mask mask indicating where the subject is
/// @param lights the light directions corresponding to each image
/// @returns an image where each pixel of the subject has been colored based on the depth
/// @pre 'images' has as many entries as 'lights' has rows
/// @pre 'lights' has 3 columns
/// @pre 'mask' and each entry in 'images' are all the same size
///
PImage computeDepths(PImage[] images, PImage mask, DenseMatrix lights) {
  int variables = mask.width * mask.height;
  int constraints = 2 * variables;
  FlexCompColMatrix A = new FlexCompColMatrix(constraints, variables);
  DenseVector b = new DenseVector(constraints);
  for(int y = 0; y < mask.height; y++) {
    for(int x = 0; x < mask.width; x++) {
      if(red(mask.get(x, y)) > 0) {
        //
        // TODO: Fill in A and b
        // We need to find V1 and V2
       // 0 = N.V1 and 0=N.V2 
        DenseVector normal = computeOneNormal(images, lights, x, y);
        
        double Nx = normal.get(X);
        double Ny = normal.get(Y);
        double Nz = normal.get(Z);
        int i = y * mask.width + x;
        
        A.set(2*i, i, Nz);
        A.set(2*i, i+1, -Nz);
        b.set(2*i, Nx);
        
        A.set(2*i+1, i, Nz);
        A.set(2*i+1, i+mask.width, -Nz);
        b.set(2*i+1, Ny);
      
        
      }
    }
  }

  //
  // Solve system
  //
  println("Computing ATA");
  int [] positions = new int[]{-mask.width, -1, 0, 1, mask.width};
  CompDiagMatrix ATA = new CompDiagMatrix(variables, variables, positions);
  for(int row = 0; row < variables; row ++) {
    for(int index = 2; index < positions.length; index++) {
      int column = row + positions[index];
      if(column >= 0 & column < variables) {
        no.uib.cipr.matrix.Vector rowVector = A.getColumn(column);
        no.uib.cipr.matrix.Vector colVector = A.getColumn(row);
        double value =  rowVector.dot(colVector);
        if(!Double.isNaN(value)) {
          ATA.set(row, column, value);
          ATA.set(column, row, value);
        }
      }
    }
  }

  println("Computing ATb");
  DenseVector ATb = new DenseVector(variables);
  A.transMult(b, ATb);

  println("Solving ATA = ATb z");
  DenseVector z = new DenseVector(variables);
  try{
    z.set(ATb);
    CG conjugateGradient = new CG(z);
    conjugateGradient.solve(ATA,ATb,z);
  } catch(IterativeSolverNotConvergedException e) {
    println("Could not solve system because of " + e.getReason());
    // Also print a tip
    println("Be sure that you have no NaN values in your matrix and vector.");
    println("Filter normals with NaN values. These can be found with Double.isNaN(.);");
    return null;
  }

  //
  // Visualize depths
  //

  // Find minimum and maximum
  double minD = Double.MAX_VALUE;
  double maxD = Double.MIN_VALUE;
  for(int y = 0; y < mask.height; y++){
    for(int x = 0; x < mask.width; x++) {
      if(red(mask.get(x,y)) > 0){
        int pixel = x + y * mask.width;
        double value = z.get(pixel);
        minD = Math.min(value, minD);
        maxD = Math.max(value, maxD);
      }
    }
  }
  // Render scaled values
  PImage result = createImage(mask.width, mask.height, RGB);
  for(int y = 0; y < mask.height; y++){
    for(int x = 0; x < mask.width; x++) {
      if(red(mask.get(x,y)) > 0) {
        int pixel = x + y * mask.width;
        double value = z.get(pixel);
        double i = (value - minD)/(maxD - minD);
        result.set(x, y, color(Math.round(255 * i)));
      }
    }
  }
  return result;
}


///
/// Loads all images that fit the given printf 'pattern' in the given 'directory'
///
/// @param directory directory to look in
/// @param pattern pattern with %d (or modification) in it somewhere
/// @returns an array of PImages for the files that met the pattern
/// @pre the 'directory' must not contain more than MAX_IMAGES images that meet the pattern
/// @pre the images must be ordered sequentially starting with 0
/// @see http://java.sun.com/j2se/1.5.0/docs/api/java/util/Formatter.html
///
PImage[] loadImageDirectory(String directory, String pattern) {
  println("Loading all images like '" + directory + SLASH + pattern + "'");
  ArrayList images = new ArrayList();
  String sanityCheck = "";
  for(int index = 0; index < MAX_IMAGES; index++){
    String filename =  String.format(pattern, index).toString();
    if(filename == sanityCheck) {
      println("Your pattern is bad! It must contain %d or some other integer placeholder!");
      return null;
    }
    PImage image = loadImage(directory + SLASH + filename);
    if(null == image) {
      if(index == 0) {
        println("No images were found, check to make sure there are images named correctly.");
      } else {
        println("Loaded " + index + " images (don't worry about the above error).");
      }
      return (PImage[])images.toArray(new PImage[0]);
    }
    images.add(image);
  }
  println("Loaded " + MAX_IMAGES + " images.");
  return (PImage[])images.toArray(new PImage[0]);
}


///
/// Reads points into a an Nx2 matrix from a file with N lines
/// File may have commas, or not, should only have 2 values (X and Y) per line
///
/// @param filename - Text file with X and Y pairs on each line
/// @pre file exists
/// @returns an Nx2 matrix with the values filled in
///
Matrix readPoints(String filename) {
  String[] lines = loadStrings(filename);
  DenseMatrix matrix = new DenseMatrix(lines.length, 2);
  for(int l = 0; l < lines.length; l++) {
    String[] values = splitTokens(lines[l], ", ");
    matrix.set(l,X,float(values[X]));
    matrix.set(l,Y,float(values[Y]));
  }
  return matrix;
}


///
/// Wrapper that displays an image
///
/// @param i image to show
///
void show(PImage i) {
  show(i, new DenseMatrix(0,0));
}


///
/// Displays an image with highlighted points
///
/// @param i image to show
/// @param points points to circle in red
///
void show(PImage i, Matrix points) {
  imageToShow = i;
  pointsToShow = points;
  
  redraw();
}


///
/// Makes small circles on the pallet image to show where the points in 'points' are
///
/// @param points points to make circles at
/// @param x horizontal offset (used for centering images, see show)
/// @param y vertical offset (used for centering images, see show)
///
void showPoints(Matrix points){
  stroke(0,0,0);
  fill(255,0,0);
  for(int r = 0; r < points.numRows(); r++) {
    ellipse((float) points.get(r,X), (float)points.get(r,Y), DOT_SIZE, DOT_SIZE);
  }
}

///
/// Reads points into a an Nx3 matrix from a file with N lines
/// File may have commas, or not, should only have 3 values (X, Y and Z) per line
///
/// @param filename - Text file with X, Y, Z tuples on each line
/// @pre file exists
/// @returns an Nx3 matrix with the values filled in
///
DenseMatrix readPoints3D(String filename) {
  String[] lines = loadStrings(filename);
  DenseMatrix matrix = new DenseMatrix(lines.length, 3);
  for(int l = 0; l < lines.length; l++) {
    String[] values = splitTokens(lines[l], "\t, ");
    matrix.set(l,X,float(values[X]));
    matrix.set(l,Y,float(values[Y]));
    matrix.set(l,Z,float(values[Z]));
  }
  return matrix;
}
