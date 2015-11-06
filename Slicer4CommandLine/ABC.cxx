

// Atlas based segmentation
// ABC outdir --inputImage img1 ... --inputImage imgn --atlasMRB atlas.mrb -b 4 --atlasFluidIters 5 --outputLabel seg

#include "ABCCLP.h"

#include "mu.h"
#include "muFile.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkMultiThreader.h"

#include "ModuleDescriptionParser.h"
#include "ModuleDescription.h"
#include "ModuleParameterGroup.h"
#include "ModuleParameter.h"

#include "vtkMRMLApplicationLogic.h"
#include "vtkMRMLNode.h"
#include "vtkMRMLScene.h"
#include "vtkSlicerApplicationLogic.h"

#include "vtkMRMLScalarVolumeNode.h"
#include "vtkMRMLSceneViewNode.h"
#include "vtkMRMLSceneViewStorageNode.h"

#include "vtkMRMLSubjectHierarchyNode.h"

#include <vtkNew.h>
#include <vtksys/SystemTools.hxx>

#include "DynArray.h"
#include "Timer.h"

// Use manually instantiated classes for the big program chunks
#define MU_MANUAL_INSTANTIATION
#include "EMSegmentationFilter.h"
#include "PairRegistrationMethod.h"
#undef MU_MANUAL_INSTANTIATION

#include "AtlasMRBRegistrationMethod.h"

#include <cstdlib>
#include <exception>
#include <iostream>

#include <string>


typedef std::vector<std::string> StringList;

int run_ABC(int argc, char** argv)
{

  PARSE_ARGS;

  if (labelImage.size() == 0 && outputImage1.size() == 0)
  {
    std::cerr << "No output label or bias-corrected image specified" << std::endl;
    return -1;
  }

  std::cout << "Locating temporary directory..." << std::endl;

  char* envTempDir = getenv("SLICER_TEMPORARY_DIR");

  if (envTempDir == 0)
    envTempDir = getenv("TEMP");
  if (envTempDir == 0)
    envTempDir = getenv("TMP");
  if (envTempDir == 0)
    envTempDir = getenv("TMPDIR");
  if (envTempDir == 0)
    envTempDir = "/tmp/Slicer";

  std::string tempDir = std::string(envTempDir);
  tempDir += std::string("/ABC-CLI-MRB");

  std::cout << "Using temporary directory: " << tempDir << std::endl;

  vtkNew<vtkSlicerApplicationLogic> appLogic;

  vtkNew<vtkMRMLScene> scene;

  vtkNew<vtkMRMLHierarchyNode> hNode;
  scene->RegisterNodeClass(hNode.GetPointer(), "vtkMRMLHierarchyNode");
  vtkNew<vtkMRMLSubjectHierarchyNode> shNode;
  scene->RegisterNodeClass(shNode.GetPointer(), "vtkMRMLSubjectHierarchyNode");

  vtkNew<vtkMRMLSceneViewNode> svNode;
  scene->RegisterNodeClass(svNode.GetPointer(), "vtkMRMLSceneViewNode");
  vtkNew<vtkMRMLSceneViewStorageNode> svsNode;
  scene->RegisterNodeClass(svsNode.GetPointer(), "vtkMRMLSceneViewStorageNode");

  appLogic->SetMRMLScene(scene.GetPointer());

  vtksys::SystemTools::MakeDirectory(tempDir);

  std::cout << "Loading MRB: " << atlasMRB << std::endl;

  bool atlasLoaded = appLogic->OpenSlicerDataBundle(
    atlasMRB.c_str(), tempDir.c_str() );

  scene->InitTraversal();

  vtkMRMLSubjectHierarchyNode* rootNode =
    dynamic_cast<vtkMRMLSubjectHierarchyNode*>(
      scene->GetNextNodeByClass("vtkMRMLSubjectHierarchyNode") );

  if (rootNode == 0)
  {
    std::cerr << "No subject hierarchy detected in MRB" << std::endl;
    vtksys::SystemTools::RemoveADirectory(tempDir);
    return -1;
  }

/*
  vtkNew<vtkCollection> childNodes;
  rootNode->GetAssociatedChildrenNodes(childNodes.GetPointer());

  vtkMRMLScalarVolumeNode* templateVol = 0;

  for (int i = 0; i < childNodes->GetNumberOfItems(); i++)
  {
    vtkMRMLScalarVolumeNode* vn = dynamic_cast<vtkMRMLScalarVolumeNode*>(
      childNodes->GetItemAsObject(i) );

    if (vn == 0)
      continue;

    if (strcmp(vn->GetName(), "template") == 0)
    {
      templateVol = vn;
      break;
    }
  }

  if (templateVol == 0)
  {
    std::cerr << "No template volume node detected" << std::endl;
    vtksys::SystemTools::RemoveADirectory(tempDir);
    return -1;
  }

  std::vector<vtkMRMLScalarVolumeNode*> priorVols;
  for (int i = 0; i < childNodes->GetNumberOfItems()-1; i++)
  {
    char s[256];
    snprintf(s, 256, "%d", i+1);

    vtkMRMLScalarVolumeNode* p_i = 0;

    for (int j = 0; j < childNodes->GetNumberOfItems(); j++)
    {
      vtkMRMLScalarVolumeNode* vn = dynamic_cast<vtkMRMLScalarVolumeNode*>(
        childNodes->GetItemAsObject(j) );

      if (vn == 0)
        continue;

      if (strcmp(vn->GetName(), s) == 0)
      {
        p_i = vn;
        break;
      }
    }

    priorVols.push_back(p_i);
  }

  if (priorVols.size() == 0)
  {
    std::cerr << "No prior volume nodes detected" << std::endl;
    return -1;
  }

std::cout << "template stor " << templateVol->GetStorageNode() << std::endl;

  std::cout << "template " << templateVol->GetName() << " ID = " << templateVol->GetID() << " file = " << templateVol->GetStorageNode()->GetFileName() << std::endl;
  for (int i = 0; i < priorVols.size(); i++)
    std::cout << "prior " << i << " " << priorVols[i]->GetName() << " ID = " << priorVols[i]->GetID() << " file = " << priorVols[i]->GetStorageNode()->GetFileName() << std::endl;
*/

  itk::MultiThreader::SetGlobalDefaultNumberOfThreads( NumberOfThreads );

  Timer* timer = new Timer();

  typedef itk::Image<unsigned char, 3> ByteImageType;
  typedef itk::Image<float, 3> FloatImageType;
  typedef itk::Image<short, 3> ShortImageType;

  typedef itk::ImageFileReader<FloatImageType> ReaderType;

  // Set up suffix string for images
  std::string suffstr = "_seg.mha";

  std::cout << "Reading input images: " << std::endl;

  DynArray<std::string> inputFiles;
  if (inputImage1.size() != 0)
    inputFiles.Append(inputImage1);
  if (inputImage2.size() != 0)
    inputFiles.Append(inputImage2);
  if (inputImage3.size() != 0)
    inputFiles.Append(inputImage3);
  if (inputImage4.size() != 0)
    inputFiles.Append(inputImage4);
  if (inputImage5.size() != 0)
    inputFiles.Append(inputImage5);

  DynArray<std::string> inputOrients;
/*
  if (inputOrient1.size() != 0)
    inputOrients.Append(inputOrient1);
  if (inputOrient2.size() != 0)
    inputOrients.Append(inputOrient2);
  if (inputOrient3.size() != 0)
    inputOrients.Append(inputOrient3);
  if (inputOrient4.size() != 0)
    inputOrients.Append(inputOrient4);
*/
  for (unsigned int i = 0; i < inputFiles.GetSize(); i++)
    inputOrients.Append(std::string("file"));

  muLogMacro(<< "Registering images using linear transform...\n");

  ByteImageType::Pointer fovmask;

  FloatImageType::Pointer templateImage;

  DynArray<FloatImageType::Pointer> images;
  DynArray<FloatImageType::Pointer> priors;

  try
  {
    typedef AtlasMRBRegistrationMethod<float, float> AtlasRegType;
    AtlasRegType::Pointer atlasreg = AtlasRegType::New();

    atlasreg->SetPrefilteringMethod(FilterMethod.c_str());
    atlasreg->SetPrefilteringIterations(FilterIterations);
    atlasreg->SetPrefilteringTimeStep(FilterTimeSteps);

    //TODO: allow parameter to set output names?
    atlasreg->SetSuffix("");

    atlasreg->SetImageFileNames(inputFiles);
    atlasreg->SetImageOrientations(inputOrients);
    // NOTE: trafo write should be disabled for Slicer
    atlasreg->SetOutputDirectory(std::string(""));

    if (atlasMapType.compare("identity") == 0)
      atlasreg->SetAtlasLinearTransformChoice(AtlasRegType::ID_TRANSFORM);
    if (atlasMapType.compare("rigid") == 0)
      atlasreg->SetAtlasLinearTransformChoice(AtlasRegType::RIGID_TRANSFORM);

    if (coregMapType.compare("identity") == 0)
      atlasreg->SetImageLinearTransformChoice(AtlasRegType::ID_TRANSFORM);
    if (coregMapType.compare("rigid") == 0)
      atlasreg->SetImageLinearTransformChoice(AtlasRegType::RIGID_TRANSFORM);

    if (RegistrationMode.compare("Fine") == 0)
      atlasreg->FastRegistrationOff();

    atlasreg->SetAtlasHierarchy(rootNode);

    // NOTE: always start from scratch in Slicer
    //muLogMacro(<< "Attempting to read previous registration results..."
    //  << std::endl);
    //atlasreg->ReadParameters();

    muLogMacro(<< "Registering and resampling images..." << std::endl);
    atlasreg->Update();

    // NOTE: Disable write, unless you have output transform node?
    //atlasreg->WriteParameters();

    fovmask = atlasreg->GetFOVMask();

    images = atlasreg->GetImages();
    priors = atlasreg->GetProbabilities();

    templateImage = atlasreg->GetAffineTemplate();
  } // end atlas reg block
  catch (...)
  {
    std::cerr << "Error while registering atlas" << std::endl;
    vtksys::SystemTools::RemoveADirectory(tempDir);
    return -1;
  }

  // Remove unpacked atlas MRB when done
  vtksys::SystemTools::RemoveADirectory(tempDir);

  // Rescale intensity of input images
  for (unsigned int k = 0; k < images.GetSize(); k++)
  {
    typedef itk::RescaleIntensityImageFilter<FloatImageType, FloatImageType>
      RescalerType;
    RescalerType::Pointer resf = RescalerType::New();
    resf->SetInput(images[k]);
    resf->SetOutputMinimum(1);
    resf->SetOutputMaximum(32000);
    resf->Update();
    images[k] = resf->GetOutput();
  }

  std:: cout << "Start segmentation..." << std::endl;
  typedef EMSegmentationFilter<FloatImageType, FloatImageType> SegFilterType;
  SegFilterType::Pointer segfilter = SegFilterType::New();

  segfilter->SetTemplateImage(templateImage);

  segfilter->SetInputImages(images);
  segfilter->SetPriors(priors);

  SegFilterType::VectorType priorweights(priorAdjustVec.size());
  for (unsigned int i = 0; i < priorAdjustVec.size(); i++)
    priorweights[i] = priorAdjustVec[i];
  segfilter->SetPriorWeights(priorweights);

  segfilter->SetInitialDistributionEstimator(InitialDistributionEstimator);

  segfilter->SetMaxBiasDegree(biasDegree);

  bool dowarp = (atlasFluidIters != 0);

  if (dowarp)
    segfilter->WarpingOn();
  else
    segfilter->WarpingOff();

  segfilter->SetWarpFluidIterations(atlasFluidIters);
  segfilter->SetWarpFluidMaxStep(atlasFluidMaxStep);

  segfilter->Update();

  std::cout << "Writing segmentation images..." << std::endl;

  {
    typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
    ShortWriterType::Pointer writer = ShortWriterType::New();

    writer->SetFileName(labelImage.c_str());
    writer->SetInput(segfilter->GetOutput());
    writer->Update();
  }

  // Write registered - bias corrected images

  DynArray<std::string> outputFiles;
  outputFiles.Append(outputImage1);
  outputFiles.Append(outputImage2);
  outputFiles.Append(outputImage3);
  outputFiles.Append(outputImage4);
  outputFiles.Append(outputImage5);

  DynArray<FloatImageType::Pointer> corrImages = segfilter->GetCorrected();

  for (unsigned int i = 0; i < corrImages.GetSize(); i++)
  {
    if (outputFiles[i].size() == 0)
      continue;

    typedef itk::RescaleIntensityImageFilter<FloatImageType, ShortImageType>
      RescalerType;
    RescalerType::Pointer resf = RescalerType::New();
    resf->SetInput(corrImages[i]);
    resf->SetOutputMinimum(1);
    resf->SetOutputMaximum(32000);
    resf->Update();

    typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
    ShortWriterType::Pointer writer = ShortWriterType::New();

    writer->SetInput(resf->GetOutput());
    writer->SetFileName(outputFiles[i].c_str());
    writer->Update();
  }
  

/*
  // Write brain class posteriors
  DynArray<ShortImageType::Pointer> posteriors =
    segfilter->GetShortPosteriors();

  for (unsigned int i = 0; i < (posteriors.GetSize()-3); i++)
  {
    std::ostringstream oss;
    oss << outputDir << mu::get_name(inputFiles[0].c_str())
      << "_posterior" << i << suffstr << std::ends;

    typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
    ShortWriterType::Pointer writer = ShortWriterType::New();

    writer->SetInput(posteriors[i]);
    writer->SetFileName(oss.str().c_str());
    writer->Update();
  }
*/

  timer->Stop();

  muLogMacro(<< "All segmentation processes took " << timer->GetElapsedHours() << " hours, ");
  muLogMacro(<< timer->GetElapsedMinutes() << " minutes, ");
  muLogMacro(<< timer->GetElapsedSeconds() << " seconds\n");

  return 0;

}

int
main(int argc, char** argv)
{

  PARSE_ARGS;

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  try
  {
    run_ABC(argc, argv);
  }
  catch (itk::ExceptionObject& e)
  {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch (std::string& s)
  {
    std::cerr << "Exception: " << s << std::endl;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cerr << "Unknown exception" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}
