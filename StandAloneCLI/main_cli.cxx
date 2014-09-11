
#include "mu.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "EMSParameters.h"
#include "EMSParametersXMLFile.h"
#include "runEMS.h"

#include <exception>
#include <iostream>

#include <getopt.h>
#include <unistd.h>

std::string minc_timestamp(int argc,char **argv)
{
  std::string timestamp;

  char cur_time[200];
  time_t t;
  struct tm *tmp;

  t = time(NULL);
  tmp = localtime(&t);

  strftime(cur_time, sizeof(cur_time), "%a %b %d %T %Y>>>", tmp);
  /* Get the time, overwriting newline */
  timestamp=cur_time;

  /* Copy the program name and arguments */
  for (int i=0; i<argc; i++) {
    timestamp+=argv[i];
    timestamp+=" ";
  }
  timestamp+="\n";

  return timestamp;
}



void
printUsage(char* progname)
{
  std::cerr << "Usage: " << progname << " <input1> [input2] [input3] <output dir> <output suffix> [options]" << std::endl
            << "Available options:" << std::endl
            << "--debug:\tdisplay debug messages" << std::endl
            << "--write-less:\tdon't write posteriors and filtered, bias corrected images"
            << std::endl
            << "--template <directory with template> - REQUIRED!"<< std::endl 
            << "--affine\t perform affine registration to the template first"<<std::endl
            << "--nl\t perform nl registration to the template "<<std::endl
            << "--filter-steps <n> use curvature-flow filter for pre-filtering images with this many steps"<<std::endl
            << "--filter-timestep <f> set time-step for curvature-flow filter"<<std::endl
            << "TODO: add more options!"<<std::endl;
}



int
main(int argc, char** argv)
{
  std::string history= minc_timestamp(argc,argv);
  int clobber=0;
  int debugflag = 0;
  int writeflag = 1;
  int affine = 0;
  int rigid_body = 0;
  int nl = 0;
  int filter_steps=1;
  double filter_timestep=0.1;

  std::string template_dir;
  
  static struct option long_options[] = {
    {"debug",       no_argument,       &debugflag, 1},
    {"writeless",   no_argument,       &writeflag, 0},
    {"template",    required_argument, 0, 't'},
    {"filter-steps",required_argument, 0, 'f'},
    {"filter-timestep",  required_argument, 0, 'F'},
    {"affine",    no_argument,       &affine, 1},
    {"rigid",     no_argument,       &rigid_body, 1},
    {"nl",        no_argument,       &nl, 1},
    {0, 0, 0, 0}
    };
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int
    option_index = 0;

    c = getopt_long (argc, argv, "t:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
    case 0:
      break;
    case 't':
      template_dir = optarg;
      break;
    case 'f':
      filter_steps=atoi(optarg);
      break;
    case 'F':
      filter_timestep=atof(optarg);
      break;
    case '?':
      /* getopt_long already printed an error message. */
    default:
      printUsage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<3 || template_dir.empty() )
  {
    printUsage(argv[0]);
    return 1;
  }

  std::string output_suffix=argv[argc-1];
  std::string output_dir=argv[argc-2];
  argc-=2;


  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  try
  {
    EMSParameters::Pointer emsp = EMSParameters::New();
    emsp->SetAtlasDirectory(template_dir);
    emsp->SetAtlasOrientation("RAI");
    
    emsp->SetOutputFormat("MINC");
    
    for(int i=optind;i<argc;i++)
      emsp->AddImage(argv[i],"RAI");

    emsp->SetFilterIterations(filter_steps);
    emsp->SetFilterTimeStep(filter_timestep);
    //emsp->SetFilterMethod();
    emsp->SetMaxBiasDegree(4);
    emsp->SetDoAtlasWarp(nl); 
    //SetAtlasWarpFluidIterations
    //SetAtlasWarpFluidMaxStep
    //SetAtlasWarpKernelWidth
    if(affine)
    {
      emsp->SetImageLinearMapType("affine");
      emsp->SetAtlasLinearMapType("affine");
    }
    else if(rigid_body)
    {
      emsp->SetImageLinearMapType("rigid");
      emsp->SetAtlasLinearMapType("rigid");
    }
    else
    {
      emsp->SetImageLinearMapType("id");
      emsp->SetAtlasLinearMapType("id");
    }
    
    for(int i=0;i<4;i++)
      emsp->AppendPriorWeight(1.0);

    emsp->SetNumberOfThreads(1);

    emsp->SetSuffix(output_suffix);
    emsp->SetOutputDirectory(output_dir);
    
    if(debugflag)
      emsp->PrintSelf(std::cerr);
    
    runEMS(emsp, debugflag, writeflag);
  }
  catch (itk::ExceptionObject& e)
  {
    std::cerr << e << std::endl;
    return -1;
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << std::endl;
    return -1;
  }
  catch (std::string& s)
  {
    std::cerr << "Exception: " << s << std::endl;
    return -1;
  }
  catch (...)
  {
    std::cerr << "Unknown exception" << std::endl;
    return -1;
  }

  return 0;

}
