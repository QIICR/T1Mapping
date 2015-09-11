#include "itkImageFileWriter.h"
#include "itkImageDuplicator.h"
#include "itkMetaDataObject.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkArray.h"


#include "itkPluginUtilities.h"
//#include "lmcurve.h"

#include "T1MappingCLP.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#define PI 3.1415926535897932384626433832795
float TR;
//double f(double t, const double *p){
//  return p[1]*exp(p[0]*t);
//}

//int fit_exp(double *par, int m_dat, double *t, double *y);

#define SimpleAttributeGetMethodMacro(name, key, type)     \
type Get##name(itk::MetaDataDictionary& dictionary)           \
{\
  type value = type(); \
  if (dictionary.HasKey(key))\
    {\
    /* attributes stored as strings */ \
    std::string valueString; \
    itk::ExposeMetaData(dictionary, key, valueString);  \
    std::stringstream valueStream(valueString); \
    valueStream >> value; \
    }\
  else\
    {\
    itkGenericExceptionMacro("Missing attribute '" key "'.");\
    }\
  return value;\
}

//SimpleAttributeGetMethodMacro(EchoTime, "MultiVolume.DICOM.EchoTime", float);
SimpleAttributeGetMethodMacro(RepetitionTime, "MultiVolume.DICOM.RepetitionTime",float);
SimpleAttributeGetMethodMacro(FlipAngle, "MultiVolume.DICOM.FlipAngle", float);

std::vector<float> GetFA(itk::MetaDataDictionary& dictionary)
{
  std::vector<float> FA;

  if (dictionary.HasKey("MultiVolume.FrameIdentifyingDICOMTagName"))
    {
    std::string tag;
    itk::ExposeMetaData(dictionary, "MultiVolume.FrameIdentifyingDICOMTagName", tag);
    if (dictionary.HasKey("MultiVolume.FrameLabels"))
      {
      // Acquisition parameters stored as text, FrameLabels are comma separated
      std::string frameLabelsString;
      itk::ExposeMetaData(dictionary, "MultiVolume.FrameLabels", frameLabelsString);
      std::stringstream frameLabelsStream(frameLabelsString);
      if (tag == "FlipAngle")
        {
        float t;
        while (frameLabelsStream >> t)
          {
          FA.push_back(t);
          frameLabelsStream.ignore(1); // skip the comma
          }
        }
      else
        {
        itkGenericExceptionMacro("Unrecognized frame identifying DICOM tag name " << tag);
        }
      }
    else
      {
      itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameLabels'.")
      }
    }
  else
    {
    itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameIdentifyingDICOMTagName'.");
    }
  
  return FA;
}


class ExpDecayCostFunction: public itk::MultipleValuedCostFunction
{
public:
  typedef ExpDecayCostFunction                    Self;
  typedef itk::MultipleValuedCostFunction   Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;
  itkNewMacro( Self );
        
  enum { SpaceDimension =  3 };
  unsigned int RangeDimension; 

  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef Superclass::MeasureType                 MeasureType, ArrayType;
  typedef Superclass::ParametersValueType         ValueType;
		      
        
  ExpDecayCostFunction()
  {
  }
        
  void SetY (const float* y, int sz) //Self signal Y
  {    
    Y.set_size (sz);
    for (int i = 0; i < sz; ++i)
      Y[i] = y[i];
    //std::cout << "Cv: " << Y << std::endl;
  }
        
  void SetX (const float* x, int sz) //Self signal X
  {
    X.set_size (sz);
    for( int i = 0; i < sz; ++i )
      X[i] = x[i];
    //std::cout << "Time: " << X << std::endl;
  }

  ArrayType GetX(){
    return X;
  }

  ArrayType GetY(){
    return Y;
  }
        
  MeasureType GetFittedValue( const ParametersType & parameters) const
  {
    MeasureType measure(RangeDimension);
    for(int i=0;i<measure.size();i++)
      {
      measure[i] = exp(-1.*X[i]*parameters[0])*parameters[1];
      }
    return measure;
  }

  MeasureType GetValue( const ParametersType & parameters) const
  {
    MeasureType measure(RangeDimension);

    for(int i=0;i<measure.size();i++)
      {
      measure[i] = Y[i]-exp(-1.*X[i]*parameters[0])*parameters[1];
      }
           
    return measure; 
  }
        
  //Not going to be used
  void GetDerivative( const ParametersType & /* parameters*/,
                      DerivativeType  & /*derivative*/ ) const
  {   
  }
        
  unsigned int GetNumberOfParameters(void) const
  {
    return 2;
  }
       
  void SetNumberOfValues(unsigned int nValues)
    {
    RangeDimension = nValues;
    }

  unsigned int GetNumberOfValues(void) const
  {
    return RangeDimension;
  }
        
protected:
  virtual ~ExpDecayCostFunction(){}
private:
        
  ArrayType X, Y;
        
  ArrayType Exponential(ArrayType X) const
  {
    ArrayType Z;
    Z.set_size(X.size());
    for (unsigned int i=0; i<X.size(); i++)
      {
      Z[i] = exp(X(i));
      }
    return Z;
  };
        
  int constraintFunc(ValueType x) const
  {
    if (x<0||x>1)
      return 1;
    else
      return 0;
            
  };
        
        
};

class DecayCostFunction: public itk::MultipleValuedCostFunction
{
public:
  typedef DecayCostFunction                    Self;
  typedef itk::MultipleValuedCostFunction   Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;
  itkNewMacro( Self );
        
  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef Superclass::MeasureType                 MeasureType, ArrayType;
  typedef Superclass::ParametersValueType         ValueType;

  enum Model {
    VFA = 0,
  };

  enum { SpaceDimension =  3 };

  DecayCostFunction()
  {
    modelType = VFA;
  }

  void SetModelType(Model mt){
    modelType = mt;
    switch(modelType){
    
    case VFA:
      initialValue = ParametersType(2);
      initialValue[0] = 0;
      //initialValue[1] = 0.0015;
      initialValue[1] = 0.02;

      // initialize parameter meaning (store this in NRRD? save units?)
      parametersMeaning.clear();
      parametersMeaning.push_back("Scale");
      parametersMeaning.push_back("T1");

      break;
    default:
      abort(); // not implemented!
    }
  }

 ParametersType GetInitialValue(){
   return initialValue;
 }

  void SetInitialValues(ParametersType initialParameters){
    // TODO: add model-specific checks
    if(initialParameters.size() != initialValue.size())
      return;
    for(int i=0;i<initialValue.size();i++)
      initialValue[i] = initialParameters[i];
  }

  void SetY (const float* y, int sz) //Self signal Y
  {    
    Y.set_size (sz);
    for (int i = 0; i < sz; ++i)
      Y[i] = y[i];
    //std::cout << "Cv: " << Y << std::endl;
  }
        
  void SetX (const float* x, int sz) //Self signal X
  {
    X.set_size (sz);
    for( int i = 0; i < sz; ++i )
      X[i] = x[i];
    //std::cout << "Time: " << X << std::endl;
  }

  ArrayType GetX(){
    return X;
  }

  ArrayType GetY(){
    return Y;
  }
        
  MeasureType GetFittedVector( const ParametersType & parameters) const
  {
    MeasureType measure(RangeDimension);

    switch(modelType){
    
    case VFA:{
      float scale = parameters[0],
          t1 = parameters[1];

      for(int i=0;i<measure.size();i++)
        {
        measure[i] =
          scale*(sin (X[i]*PI/180)*(1-exp(-1.*(TR/t1))))/(1-(exp(-1.*(TR/t1)))*cos (X[i]*PI/180));
          
        }
      break;
    }
    default:
      abort();
    }

    return measure;
  }

  float GetFittedValue(const ParametersType & parameters, float x) const
  {
    float measure;
    switch(modelType){
    
    case VFA:{
      float scale = parameters[0],
          t1 = parameters[1];
      measure =
        scale*(sin (x*PI/180)*(1-exp(-1.*(TR/t1))))/(1-(exp(-1.*(TR/t1)))*cos (x*PI/180));  
        
      break;
    }
    default:
      abort();
    }

    return measure;
  }

  MeasureType GetValue( const ParametersType & parameters) const
  {
    MeasureType measure(RangeDimension);

    switch(modelType){
    
    case VFA:{
      float scale = parameters[0],
          t1 = parameters[1];

      for(int i=0;i<measure.size();i++)
        {
        measure[i] =
          
          Y[i]-scale*(sin (X[i]*PI/180)*(1-exp(-1.*(TR/t1))))/(1-(exp(-1.*(TR/t1)))*cos (X[i]*PI/180));
        }
      break;
    }
    default:
      abort(); // not implemented
    }
           
    return measure; 
  }
        
  //Not going to be used
  void GetDerivative( const ParametersType & /* parameters*/,
                      DerivativeType  & /*derivative*/ ) const
  {   
  }
        
  unsigned int GetNumberOfParameters(void) const
  {
    switch(modelType){
    case VFA: return 2;
    default: return 0; // should never get here
    }
  }
       
  void SetNumberOfValues(unsigned int nValues)
    {
    RangeDimension = nValues;
    }

  unsigned int GetNumberOfValues(void) const
  {
    return RangeDimension;
  }

  Model GetModelType() const {
    return modelType;
  }
        
protected:
  virtual ~DecayCostFunction(){}
private:
        
  ArrayType X, Y;

  unsigned int RangeDimension;
  Model modelType;
  ParametersType initialValue;
  std::vector<std::string> parametersMeaning;
        
  int constraintFunc(ValueType x) const
  {
    if (x<0||x>1)
      return 1;
    else
      return 0;
            
  };
        
        
};

// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
void OnlineVariance(itk::MultipleValuedCostFunction::MeasureType &values,
    double &mean, double &SD){
  double n = 0, M2 = 0;
     
  for(unsigned int i=0;i<values.GetSize();i++){
    double x = values[i];
    n++;
    double delta = x - mean;
    mean = mean + delta/n;
    M2 = M2 + delta*(x - mean);
  }
  SD = sqrt(M2/(n-1));
}


const   unsigned int VectorVolumeDimension = 3;
typedef float                                                 VectorVolumePixelType;
typedef itk::VectorImage<VectorVolumePixelType, VectorVolumeDimension> VectorVolumeType;
typedef VectorVolumeType::RegionType              VectorVolumeRegionType;
typedef itk::ImageFileReader<VectorVolumeType>             VectorVolumeReaderType;

typedef unsigned char                                        MaskVolumePixelType;
typedef float                                                MapVolumePixelType;
typedef itk::Image<MaskVolumePixelType, 3>                   MaskVolumeType;
typedef itk::Image<MapVolumePixelType, 3>                    MapVolumeType;
typedef itk::ImageFileReader<MaskVolumeType>                 MaskVolumeReaderType;
typedef itk::ImageFileWriter<MapVolumeType>                  MapWriterType;
typedef itk::ImageFileWriter<VectorVolumeType>               FittedVolumeWriterType;

typedef itk::Image<float,VectorVolumeDimension> OutputVolumeType;
typedef itk::ImageDuplicator<VectorVolumeType> DuplicatorType;
typedef itk::ImageFileWriter< MapVolumeType> MapVolumeWriterType;

typedef itk::ImageRegionConstIterator<VectorVolumeType> InputVectorVolumeIteratorType;
typedef itk::ImageRegionIterator<VectorVolumeType> OutputVectorVolumeIteratorType;
typedef itk::ImageRegionConstIterator<MaskVolumeType> MaskVolumeIteratorType;
typedef itk::ImageRegionIterator<MapVolumeType> MapVolumeIteratorType;

void SaveMap(MapVolumeType::Pointer map, std::string fileName){
  MapWriterType::Pointer writer = MapWriterType::New();
  writer->SetInput(map);
  writer->SetFileName(fileName.c_str());
  writer->SetUseCompression(1);
  writer->Update();
}

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
int main( int argc, char * argv[])
{
  PARSE_ARGS;

  if((FAToInclude.size()>0) && (FAToExclude.size()>0)){
    std::cerr << "ERROR: Either inclusion or exclusion flip-angles list can be specified, not both!" << std::endl;
    return -1;
  }

  //Read VectorVolume
  VectorVolumeReaderType::Pointer multiVolumeReader 
    = VectorVolumeReaderType::New();
  multiVolumeReader->SetFileName(imageName.c_str() );
  multiVolumeReader->Update();
  VectorVolumeType::Pointer inputVectorVolume = multiVolumeReader->GetOutput();

  // Read mask
  MaskVolumeType::Pointer maskVolume;
  if(maskName != ""){
    MaskVolumeReaderType::Pointer maskReader = MaskVolumeReaderType::New();
    maskReader->SetFileName(maskName.c_str());
    maskReader->Update();
    maskVolume = maskReader->GetOutput();
  } else {
    maskVolume = MaskVolumeType::New();
    maskVolume->SetRegions(inputVectorVolume->GetLargestPossibleRegion());
    maskVolume->CopyInformation(inputVectorVolume);
    maskVolume->Allocate();
    maskVolume->FillBuffer(1);
  }

  //Look for tags representing the acquisition parameters
  //
  //

  // Trigger times
  std::vector<float> FA;
  // list of flip-angles to be passed to the optimizer
  float *FAPtr, *imageValuesPtr;
  // "true" for the flip-angle and measurement pair to be used in fitting
  bool *FAMask;
  int FATotal, FASelected;
  try {

    FA = GetFA(inputVectorVolume->GetMetaDataDictionary());
    FATotal = FA.size();
    FAMask = new bool[FATotal];

    // if the inclusion list is non-empty, use only the values requested by the
    // user
    if(FAToInclude.size()){
      memset(FAMask,false,sizeof(bool)*FATotal);
      FASelected = 0;
      for(int i=0;i<FATotal;i++){
        if(std::find(FAToInclude.begin(), FAToInclude.end(), FA[i]) 
            != FAToInclude.end()){
          FAMask[i] = true;
          FASelected++;
        }
      }
    // if the exclusion list is non-empty, do not use the values requested by the
    // user
    } else if(FAToExclude.size()) {
      memset(FAMask,true,sizeof(bool)*FATotal);
      FASelected = FATotal;
      for(int i=0;i<FATotal;i++){
        if(std::find(FAToExclude.begin(), FAToExclude.end(), FA[i]) 
            != FAToExclude.end()){
          FAMask[i] = false;
          FASelected--;
        }
      }
    } else {
      // by default, all flip-angles will be used
      FASelected = FATotal;
      memset(FAMask,true,sizeof(bool)*FATotal);
    }

    if(FASelected<2){
      std::cerr << "ERROR: Less than 2 values selected, cannot do the fit!" << std::endl;
      return -1;
    }

    FAPtr = new float[FASelected];
    imageValuesPtr = new float[FASelected];
    int j = 0;
    std::cout << "Will use the following flip-angles: ";
    for(int i=0;i<FATotal;i++){
      if(FAMask[i]){
        std::cout << FA[i] << " ";
        FAPtr[j++] = FA[i];
      }
    }
    std::cout << std::endl;
  } catch (itk::ExceptionObject &exc) {
    itkGenericExceptionMacro(<< exc.GetDescription() 
            << " Image " << imageName 
            << " does not contain sufficient attributes to support algorithms.");
    return EXIT_FAILURE;
  }

  // RepetitionTime
  //float TR = 0.0;
  try 
    {
    TR = GetRepetitionTime(inputVectorVolume->GetMetaDataDictionary())/1000;
    }
  catch (itk::ExceptionObject &exc)
    {
   itkGenericExceptionMacro(<< exc.GetDescription() 
            << " Image " << imageName.c_str() 
            << " does not contain sufficient attributes to support algorithms.");
    return EXIT_FAILURE;
    }

  // allocate output maps
  DecayCostFunction::Model modelType;
  std::vector<MapVolumeType::Pointer> parameterMapVector;
  std::vector<MapVolumeIteratorType> parameterMapItVector;

  if(modelName == "VFA")
    modelType = DecayCostFunction::VFA;
  else {
    std::cerr << "ERROR: Unknown model type specified!" << std::endl;
    return -1;
  }

  DecayCostFunction::Pointer costFunction = DecayCostFunction::New();
  costFunction->SetModelType(modelType);

  parameterMapVector.resize(costFunction->GetNumberOfParameters());
  for(int i=0;i<costFunction->GetNumberOfParameters();i++){
    parameterMapVector[i] = MapVolumeType::New();
    parameterMapVector[i]->SetRegions(maskVolume->GetLargestPossibleRegion());
    parameterMapVector[i]->Allocate();
    parameterMapVector[i]->FillBuffer(0);
    // note mask is initialized even if not passed by the user
    parameterMapVector[i]->CopyInformation(maskVolume);
    parameterMapVector[i]->FillBuffer(0);

    parameterMapItVector.push_back(
          MapVolumeIteratorType(parameterMapVector[i],parameterMapVector[i]->GetLargestPossibleRegion()));
    parameterMapItVector[i].GoToBegin();
  }

  // R^2 and fitted values volumes are calculated independently of the model
  MapVolumeType::Pointer rsqrMap = MapVolumeType::New();
  rsqrMap->SetRegions(maskVolume->GetLargestPossibleRegion());
  rsqrMap->Allocate();
  rsqrMap->FillBuffer(0);
  rsqrMap->CopyInformation(maskVolume);
  rsqrMap->FillBuffer(0);

  DuplicatorType::Pointer dup = DuplicatorType::New();
  dup->SetInputImage(inputVectorVolume);
  dup->Update();
  VectorVolumeType::Pointer fittedVolume = dup->GetOutput();
  VectorVolumeType::PixelType zero = VectorVolumeType::PixelType(FA.size());
  for(int i=0;i<FA.size();i++)
    zero[i] = 0;
  fittedVolume->FillBuffer(zero);

  InputVectorVolumeIteratorType vvIt(inputVectorVolume, inputVectorVolume->GetLargestPossibleRegion());
  OutputVectorVolumeIteratorType fittedIt(fittedVolume, fittedVolume->GetLargestPossibleRegion());

  MaskVolumeIteratorType mvIt(maskVolume, maskVolume->GetLargestPossibleRegion());
  MapVolumeIteratorType rsqrIt(rsqrMap, rsqrMap->GetLargestPossibleRegion());

  itk::LevenbergMarquardtOptimizer::Pointer optimizer = itk::LevenbergMarquardtOptimizer::New();

  int cnt = 0;

  for(vvIt.GoToBegin();!vvIt.IsAtEnd();++vvIt){
    //if(cnt>10)
    //  break;
    VectorVolumeType::IndexType index = vvIt.GetIndex();
    VectorVolumeType::PixelType vectorVoxel = vvIt.Get();
    VectorVolumeType::PixelType fittedVoxel(vectorVoxel.GetSize());
    for(int i=0;i<fittedVoxel.GetSize();i++)
      fittedVoxel[i] = 0;

    if(mvIt.Get() && vectorVoxel[0]){
      //cnt++;

      // use only those values that were requested by the user
      costFunction->SetX(FAPtr, FASelected);
      const float* imageVector = const_cast<float*>(vectorVoxel.GetDataPointer());
      int j = 0;
      for(int i=0;i<FATotal;i++){
        if(FAMask[i]){
          imageValuesPtr[j++] = imageVector[i];
        }
      }

      int numberOfSelectedPoints = j;
      costFunction->SetNumberOfValues(numberOfSelectedPoints);

      costFunction->SetY(imageValuesPtr,FASelected);

      DecayCostFunction::ParametersType initialValue = costFunction->GetInitialValue();
      initialValue[0] = vectorVoxel[0];
      DecayCostFunction::MeasureType temp = costFunction->GetValue(initialValue);
      optimizer->UseCostFunctionGradientOff();
      optimizer->SetCostFunction(costFunction);

      itk::LevenbergMarquardtOptimizer::InternalOptimizerType *vnlOptimizer = optimizer->GetOptimizer();
      vnlOptimizer->set_f_tolerance(1e-4f);
      vnlOptimizer->set_g_tolerance(1e-4f);
      vnlOptimizer->set_x_tolerance(1e-5f);
      vnlOptimizer->set_epsilon_function(1e-9f);
      vnlOptimizer->set_max_function_evals(200);

      try {
        optimizer->SetInitialPosition(initialValue);
        optimizer->StartOptimization();
      } catch(itk::ExceptionObject &e) {
        std::cerr << " Exception caught: " << e << std::endl;
      }

      itk::LevenbergMarquardtOptimizer::ParametersType finalPosition;

      finalPosition = optimizer->GetCurrentPosition();
      for(int i=0;i<fittedVoxel.GetSize();i++){
        fittedVoxel[i] = costFunction->GetFittedValue(finalPosition, FA[i]);
        //std::cerr << "Final position: " << finalPosition << std::endl;
        //std::cout << fittedVoxel[i] << " ";
      }

      //std::cout << std::endl;
      fittedIt.Set(fittedVoxel);
      //std::cout << "Fitted voxel: " << fittedVoxel << " params: " << finalPosition << std::endl;
      switch(modelType){
      case DecayCostFunction::VFA:{
        parameterMapItVector[0].Set(finalPosition[0]);
        parameterMapItVector[1].Set(finalPosition[1]*1e+3);
        break;
      }

      default: abort();
      }

      // initialize the rsqr map
      // see PkModeling/CLI/itkConcentrationToQuantitativeImageFilter.hxx:452
      {
        double rms = optimizer->GetOptimizer()->get_end_error();
        double SSerr = rms*rms*FASelected;
        double sumSquared = 0.0;
        double sum = 0.0;
        double rSquared = 0.0;

        for (unsigned int i=0; i < FASelected; ++i){
          sum += imageValuesPtr[i];
          sumSquared += (imageValuesPtr[i]*imageValuesPtr[i]);
        }
        double SStot = sumSquared - sum*sum/(double)FASelected;
        rSquared = 1.0 - (SSerr / SStot);
        rsqrIt.Set(rSquared);
      }
    }

    for(int i=0;i<costFunction->GetNumberOfParameters();i++){
      ++parameterMapItVector[i];
    }
    ++rsqrIt;++mvIt;++fittedIt;
  }

  switch(modelType){
  case DecayCostFunction::VFA:{
    if(T1MapFileName.size())
      SaveMap(parameterMapVector[1], T1MapFileName);
    break;
  }
  default:abort();
  }

  if(rsqrVolumeFileName.size())
    SaveMap(rsqrMap, rsqrVolumeFileName);

  if(fittedVolumeFileName.size()){
    FittedVolumeWriterType::Pointer writer = FittedVolumeWriterType::New();
    fittedVolume->SetMetaDataDictionary(inputVectorVolume->GetMetaDataDictionary());
    writer->SetInput(fittedVolume);
    writer->SetFileName(fittedVolumeFileName.c_str());
    writer->SetUseCompression(1);
    writer->Update();
  }

  return EXIT_SUCCESS;
}
