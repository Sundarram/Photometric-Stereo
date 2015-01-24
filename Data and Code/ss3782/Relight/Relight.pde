///
/// CS-583: Introduction to Computer Vision
/// https://www.cs.drexel.edu/~kon/introcompvis/
///
/// Project 3 - Photometric Stereo - Relighting Applet
///
/// @note This was intended to be an applet, but it still uses too much memory.
///       To minimize memory albedos and normals are stored for only the pixels that have a white
///       mask. Throughout the code a 'pixelCount' variable is used in loops to index into the list
///       of these values. The entire set of un-masked values are stored so that relighting is fast.
///
/// @note The items marked TODO must be completed by you.
///

// Needed for reading files
import java.io.FileNotFoundException;
import java.util.*;
// Needed for the various Matrix and Vector operations
import no.uib.cipr.matrix.*;
// Needed for functions like Max and Min
import java.lang.Math.*;
// Needed for the .ini file parsing
import org.ini4j.*;

//
// Constants
//


// Homographic indexes
static final int X = 0;
static final int Y = 1;
static final int Z = 2;

// Whether or not to print debug information
static final boolean DEBUG = true;

// System constants
static final String SLASH = System.getProperty("file.separator");
static final int MAX_IMAGES = 12;

// ini file reading constants including filename and key names
static final String DEFAULT_INI_FILE = "photometricStereo.ini";
static final String CHROME = "Chrome";
static final String PATH = "path";
static final String OUTPUT = "output";
static final String SUBJECTS = "Subjects";
static final String MASK = ".mask.png";
static final String PATTERN = ".%d.png";

//
// Global variables (computed by 'setup')
//
PImage mask;
// Bounding box of mask's white region
int minX;
int maxX;
int minY;
int maxY;
// (number of subject pixels) x 3 matrices to store all the albedos and normals
DenseMatrix albedos;
DenseMatrix normals;
boolean uniformAlbedo = false;



///
/// Set up screen size, read files, initialize globals
///
void setup() {
  size(500, 300);
  load(0);
}


///
/// Calls 'load' to reload a new section when a key (1-9) is pressed
///
void keyPressed() {
  if(' ' == key) {
    uniformAlbedo = !uniformAlbedo;
  } else {
    int section = key - '1';
    load(section);
  }
}

///
/// Loads the proper image set into the global variables
///
/// @param index which path (as indexed in the ini file) to load
///
void load(int index) {
  // Read .ini file
  Ini ini = new Ini();
  try {
    BufferedReader reader = createReader(DEFAULT_INI_FILE);
    if(reader == null) {
      throw new FileNotFoundException("user specified file");
    }
    ini.load(reader);
    String path = ini.get(SUBJECTS).get(PATH, index);
    String lightPath = ini.get(CHROME).get(OUTPUT);
    println("LOADING...");
    load(path, lightPath);
  } catch(FileNotFoundException e) {
    println("You must have the file '" + e.getMessage() + "' present.");
  } catch(IOException e) {
    println("Failure reading ini file: " + e.getMessage());
    exit();
  } catch(IndexOutOfBoundsException e) {
    println("No such image set!");
  }
}



///
/// Loads the model from the images in 'path' given the light directions in the file 'lightPath'
///
void load(String path, String lightPath) {
  //
  // Load images, light sources
  //
  PImage[] images = loadImageDirectory(path, path + PATTERN);
  
  DenseMatrix lights = readPoints3D(lightPath);
 

  mask = loadImage(path + SLASH + path + MASK);

  //
  // TODO: Compute bounding box for subject from mask by setting minX,minY,maxX,maxY from mask
  // TODO: Count number of white pixels in mask
  //
  // Note: we do this to save space, and ensure the subject is centered in the window
  //
  int pixelCount = 0;
  minX = Integer.MAX_VALUE;
  maxX = Integer.MIN_VALUE;
  minY = Integer.MAX_VALUE;
  maxY = Integer.MIN_VALUE;
  
  for (int y = 0; y < mask.height; y++) {
    for (int x = 0; x < mask.width; x++) {
      if (red(mask.get(x, y)) > 0) {
        minX = min(minX, x);
        minY = min(minY, y);
        maxX = max(maxX, x);
        maxY = max(maxY, y);
         pixelCount+=1;
      }
    }
  }



  //
  // TODO: Compute all normals and save in a [pixelCount x 3] matrix - 'normals'
  // TODO: Compute all albedos and save in a [pixelCount x 3] matrix - 'albedos'
  //
  int subjectWidth = maxX - minX + 1;
  int subjectHeight = maxY - minY + 1;
  normals = new DenseMatrix(pixelCount, 3);
  albedos = new DenseMatrix(pixelCount, 3);
  // Keep track of pixel index in mask (note, we do all this to minimize the waste of space)
  // For each pixel that has a white mask value, calculate the normal and albedo
  // Hint: use red(mask.get(x,y)) > 0  to test if the mask is white
  int whitecnt = 0;

  for (int y = minY; y <= maxY; y++) {
    for (int x = minX; x <= maxX; x++) {
      if (red(mask.get(x, y)) > 0) {
        
        DenseVector norm = computeOneNormal(images, lights, x, y);
        DenseVector alb = computeOneAlbedo(images, lights, x, y, norm);
       
        for (int j = 0; j < 3; j++) {
          normals.set(whitecnt, j, norm.get(j));
          albedos.set(whitecnt, j, alb.get(j));
        }
        whitecnt+=1;
      }
    }
  }

}



//
// TODO: Copy 'computeOneAlbedo' and 'computeOneNormal' from your successful Stereo.pde
//       since they are used above
//
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
/// Main function
///
void draw() {
  background(0);
  // Make up a light diroction from mouse
  float theta = map(mouseX, 0, width, -PI/2, PI/2);
  float phi = map(mouseY, 0, height, -PI/2, PI/2);
 
  DenseVector lightDirection = new DenseVector( new double[]{sin(theta) * cos(phi), cos(theta) * sin(phi), cos(theta)});
  lightDirection.scale(1/lightDirection.norm(DenseVector.Norm.Two));
  int pixelCount = 0;
  for(int y = minY; y <= maxY; y++) {
    for(int x = minX; x <= maxX; x++) {
      if(red(mask.get(x,y)) == 0) {
        continue;
      }
      //
      // TODO: Compute the color at (x,y)
      //
      //color c;
      //brightness = dot product of surfacenormal and light direction
     DenseVector brightness;

      DenseVector N= new DenseVector(3);
      N.set(X,normals.get(pixelCount, X));
      N.set(Y,normals.get(pixelCount, Y));
      N.set(Z,normals.get(pixelCount, Z));
    
   
      
      if (uniformAlbedo) 
      {
        brightness = new DenseVector(new double[] {0.75, 0.75, 0.75});
      } else 
      {
           brightness = new DenseVector(3);
           brightness.set(X, albedos.get(pixelCount, X));
           brightness.set(Y,albedos.get(pixelCount, Y));
            brightness.set(Z,albedos.get(pixelCount, Z));

      }
// normal and lightDirection are both unit vectors.
//so the dot product is 0-1 instead of 0-255.
   
   //   intensity.scale(normal*lightDirection*255.0);
      double dotproduct =N.dot(lightDirection);
      brightness.scale(dotproduct*255.0);
      
      color c= color((int)(brightness.get(X)), ((int)brightness.get(Y)), (int)(brightness.get(Z)));
      pixelCount+=1;



      // Center the object by shifting which pixel we set
      int i = (x - minX) + (width - (maxX - minX))/2;
      int j = (y - minY) + (height - (maxY - minY))/2;
      set(i,j, c);
    }
  }
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
      println("Your pattern is bad! It must contain %1$d or some other integer placeholder!");
      return null;
    }
    PImage image = loadImage(directory + SLASH + filename);
    if(null == image) {
      println("Loaded " + index + " images.");
      return (PImage[])images.toArray(new PImage[0]);
    }
    images.add(image);
  }
  println("Loaded " + MAX_IMAGES + " images.");
  return (PImage[])images.toArray(new PImage[0]);
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
