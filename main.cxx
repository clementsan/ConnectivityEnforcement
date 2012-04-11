#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h> 
#include <itkImageFileWriter.h>

using namespace itk;
using namespace std;

const int Dimension = 3;
typedef unsigned short ImagePixelType;
typedef Image<ImagePixelType,Dimension>  ImageType;

typedef ImageFileReader<ImageType> VolumeReaderType;
typedef ImageFileWriter<ImageType> VolumeWriterType;


static int NoDiagConnect (unsigned short *image, int *dim);
static void clear_edge(unsigned short *image, int *dims, int clear_label);


int main(int argc, const char* argv[])
{

  ImageType::Pointer image ;
  VolumeReaderType::Pointer labelReader = VolumeReaderType::New();

  // Reading image
  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  imageReader->SetFileName(argv[1]) ;
  try
    {
      imageReader->Update() ;
    }
  catch (ExceptionObject err)
    {
      cerr<<"Exception object caught!"<<std::endl;
      cerr<<err<<std::endl;
      exit(0) ;
    }
  image = imageReader->GetOutput() ;

  ImageType::IndexType nullIndex;
  nullIndex[0] = 0;
  nullIndex[1] = 0;
  nullIndex[2] = 0;
  
  ImagePixelType *data = &((*image)[nullIndex]);
  ImageType::RegionType imageRegion = image->GetBufferedRegion();
  int dim[3];
  dim[0] = imageRegion.GetSize(0);
  dim[1] = imageRegion.GetSize(1);
  dim[2] = imageRegion.GetSize(2);

  clear_edge(data, dim, 0);
  NoDiagConnect(data,dim);

  VolumeWriterType::Pointer writer = VolumeWriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(image);
  try
    {
      writer->Update();
    }
    catch (itk::ExceptionObject & err)
      {
	std::cerr<<"Exception object caught!"<<std::endl;
	std::cerr<<err<<std::endl;
      }
}


static int
NoDiagConnect (unsigned short *image, int *dim) 
  // does not allow connection via diagonals only, enforces strict 6 connectedness
  // image has to be of type unsigned short
{

  //z axis
  int dimx = dim[0];
  int dimy = dim[1];
  int dimz = dim[2];
  bool correctionNeeded = true;
  int cnt = 0;

  while (correctionNeeded) {
    cnt++;
    //if (debug) cout << "NoDiag scan " << cnt << endl; 
    correctionNeeded = false;
    int dy = dimx*dimy;
    int dx = dimx;

   for (int i = 1; i < dimx - 1; i++) {
    for (int j = 1; j < dimy - 1; j++) {
      for (int k = 1; k < dimz - 1; k++) {
       unsigned short val = image[i + j * dimx + k * dy];
       if (val != 0) {
         // x,y 
         if ((image[i-1+j*dx+k*dy] == 0) && (image[i+(j-1)*dx+k*dy] == 0) && (image[i-1+(j-1)*dx+k*dy] != 0)) {
           correctionNeeded = true;
           image[i-1+j*dx+k*dy] = val;
         }
         if ((image[i+1+j*dx+k*dy] == 0) && (image[i+(j+1)*dx+k*dy] == 0) && (image[i+1+(j+1)*dx+k*dy] != 0)) {
           correctionNeeded = true;
           image[i+1+j*dx+k*dy] = val;
         }
         if ((image[i+1+j*dx+k*dy] == 0) && (image[i+(j-1)*dx+k*dy] == 0) && (image[i+1+(j-1)*dx+k*dy] != 0)) {
           correctionNeeded = true;
           image[i+1+j*dx+k*dy] = val;
         }
         if ((image[i-1+j*dx+k*dy] == 0) && (image[i+(j+1)*dx+k*dy] == 0) && (image[i-1+(j+1)*dx+k*dy] != 0)) {
           correctionNeeded = true;
           image[i-1+j*dx+k*dy] = val;
         }
         
         // xz
         if ((image[i-1+j*dx+k*dy] == 0) && (image[i+j*dx+(k-1)*dy] == 0) && (image[i-1+j*dx+(k-1)*dy] != 0)) {
           correctionNeeded = true;
           image[i-1+j*dx+k*dy] = val;
         }
         if ((image[i+1+j*dx+k*dy] == 0) && (image[i+j*dx+(k+1)*dy] == 0) && (image[i+1+j*dx+(k+1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+1+j*dx+k*dy] = val;
         }
         if ((image[i+1+j*dx+k*dy] == 0) && (image[i+j*dx+(k-1)*dy] == 0) && (image[i+1+j*dx+(k-1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+1+j*dx+k*dy] = val;
         }
         if ((image[i-1+j*dx+k*dy] == 0) && (image[i+j*dx+(k+1)*dy] == 0) && (image[i-1+j*dx+(k+1)*dy] != 0)) {
           correctionNeeded = true;
           image[i-1+j*dx+k*dy] = val;
         }
         
         // yz
         if ((image[i+(j-1)*dx+k*dy] == 0) && (image[i+j*dx+(k-1)*dy] == 0) && (image[i+(j-1)*dx+(k-1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+(j-1)*dx+k*dy] = val;
         }
         if ((image[i+(j+1)*dx+k*dy] == 0) && (image[i+j*dx+(k+1)*dy] == 0) && (image[i+(j+1)*dx+(k+1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+(j+1)*dx+k*dy] = val;
         }
         if ((image[i+(j+1)*dx+k*dy] == 0) && (image[i+j*dx+(k-1)*dy] == 0) && (image[i+(j+1)*dx+(k-1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+(j+1)*dx+k*dy] = val;
         }
         if ((image[i+(j-1)*dx+k*dy] == 0) && (image[i+j*dx+(k+1)*dy] == 0) && (image[i+(j-1)*dx+(k+1)*dy] != 0)) {
           correctionNeeded = true;
           image[i+(j-1)*dx+k*dy] = val;
         }
       }
     }
      }
    }
  }

  return 1;
}

static void clear_edge(unsigned short *image, int *dims, int clear_label)
  // clears the edge of the image
{
  int size_plane = dims[0]*dims[1];
  int size_line = dims[0];

  for (int z = 0; z < dims[2]; z++) {
    for (int y = 0; y < dims[1]; y++) {
      if ( (y == 0) || (y == dims[1]-1) ||
        (z == 0) || (z == dims[2]-1) ) { // draw whole plane
        for (int x = 0; x < dims[0] ; x++) 
          image[x +  size_line * y + size_plane * z] = clear_label;
        } else { // draw edges of x
        image[0 +  size_line * y + size_plane * z] = clear_label;
        image[size_line - 1 +  size_line * y + size_plane * z] = clear_label;
      }
    }
  }

}
