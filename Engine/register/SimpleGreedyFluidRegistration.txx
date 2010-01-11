
#ifndef _SimpleGreedyFluidRegistration_txx
#define _SimpleGreedyFluidRegistration_txx

#include "itkDerivativeImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkWarpImageFilter.h"
#include "itkWarpVectorImageFilter.h"

#include "VectorBlurImageFilter.h"

template <class TPixel, unsigned int Dimension>
SimpleGreedyFluidRegistration<TPixel, Dimension>
::SimpleGreedyFluidRegistration()
{
  m_Iterations = 20;
  m_TimeStep = 0.05;
  m_KernelWidth = 1.0;
  m_Modified = false;
}

template <class TPixel, unsigned int Dimension>
SimpleGreedyFluidRegistration<TPixel, Dimension>
::~SimpleGreedyFluidRegistration()
{
}

template <class TPixel, unsigned int Dimension>
void
SimpleGreedyFluidRegistration<TPixel, Dimension>
::SetFixedImages(const std::vector<ImagePointer>& images)
{
  m_FixedImages = images;
  m_Modified = true;
}

template <class TPixel, unsigned int Dimension>
void
SimpleGreedyFluidRegistration<TPixel, Dimension>
::SetMovingImages(const std::vector<ImagePointer>& images)
{
  m_MovingImages = images;
  m_Modified = true;
}

template <class TPixel, unsigned int Dimension>
void
SimpleGreedyFluidRegistration<TPixel, Dimension>
::Update()
{
  if (!m_Modified)
    return;

  // Initialize out = input moving
  m_OutputImages.clear();
  for (unsigned int c = 0; c < m_MovingImages.size(); c++)
  {
    typedef itk::ImageDuplicator<ImageType> DuperType;
    typename DuperType::Pointer dupef = DuperType::New();
    dupef->SetInputImage(m_MovingImages[c]);
    dupef->Update();
    m_OutputImages.push_back(dupef->GetOutput());
  }

  // Initialize phi(x) = x
  m_DeformationField = DeformationFieldType::New();
  m_DeformationField->CopyInformation(m_FixedImages[0]);
  m_DeformationField->Allocate();

  DisplacementType zerov;
  zerov.Fill(0.0);
  m_DeformationField->FillBuffer(zerov);

  typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType> IteratorType;
  IteratorType it(m_DeformationField, m_DeformationField->GetLargestPossibleRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    ImagePointType p;
    m_DeformationField->TransformIndexToPhysicalPoint(it.GetIndex(), p);

    DisplacementType v;
    for (unsigned int i = 0; i < Dimension; i++)
      v[i] = p[i];
    it.Set(v);
  }

  for (unsigned int i = 0; i < m_Iterations; i++)
    this->Step();
}

template <class TPixel, unsigned int Dimension>
void
SimpleGreedyFluidRegistration<TPixel, Dimension>
::Step()
{
  unsigned int numChannels = m_FixedImages.size();

  DisplacementType zerov;
  zerov.Fill(0.0);

  // Compute velocity field
  // v = sum_c { (fixed_c - moving_c) * grad(moving_c) }
  DeformationFieldPointer velocF = DeformationFieldType::New();
  velocF->CopyInformation(m_FixedImages[0]);
  velocF->Allocate();
  velocF->FillBuffer(zerov);

  typedef itk::ImageRegionIteratorWithIndex<DeformationFieldType> IteratorType;
  IteratorType it(velocF, velocF->GetLargestPossibleRegion());

  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
  {
    typedef itk::DerivativeImageFilter<ImageType, ImageType>
      DerivativeFilterType;

    for (unsigned int dim = 0; dim < Dimension; dim++)
    {
      typename DerivativeFilterType::Pointer derivf = DerivativeFilterType::New();
      derivf->SetInput(m_OutputImages[ichan]);
      derivf->SetDirection(dim);
      derivf->SetUseImageSpacingOn();
      derivf->Update();

      ImagePointer grad_d = derivf->GetOutput();

      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        DisplacementType v = it.Get();

        ImageIndexType ind = it.GetIndex();
        double w =
          m_FixedImages[ichan]->GetPixel(ind) -
          m_OutputImages[ichan]->GetPixel(ind);

        v[dim] += w * grad_d->GetPixel(ind) * m_TimeStep;

        it.Set(v);
      }
    } // for dim
  } // for ichan

  // Apply Green's kernel to velocity field
  typedef VectorBlurImageFilter<DeformationFieldType, DeformationFieldType>
    DeformationSmootherType;

  typename DeformationSmootherType::Pointer defsmoother = DeformationSmootherType::New();
  typename DeformationFieldType::SizeType radii;
  radii.Fill(5);
  defsmoother->SetRadius(radii);
  defsmoother->SetVariance(m_KernelWidth*m_KernelWidth);
  defsmoother->SetInput(velocF);
  defsmoother->Update();

  velocF = defsmoother->GetOutput();

  // Compose velocity field
  // h(x) = h( g(x) ) where g(x) = x + v
  typedef itk::WarpVectorImageFilter<
    DeformationFieldType, DeformationFieldType, DeformationFieldType>
    ComposerType;
  typename ComposerType::Pointer compf = ComposerType::New();
  compf->SetInput(m_DeformationField);
  compf->SetDeformationField(velocF);
  compf->SetOutputDirection(m_DeformationField->GetDirection());
  compf->SetOutputOrigin(m_DeformationField->GetOrigin());
  compf->SetOutputSpacing(m_DeformationField->GetSpacing());
  compf->Update();

  m_DeformationField = compf->GetOutput();

  DeformationFieldPointer dispF = defsmoother->GetOutput();
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    ImageIndexType ind = it.GetIndex();

    ImagePointType p;
    m_DeformationField->TransformIndexToPhysicalPoint(ind, p);

    DisplacementType u;
    for (unsigned int i = 0; i < Dimension; i++)
      u[i] = p[i];
    dispF->SetPixel(ind, m_DeformationField->GetPixel(ind) - u);
  }

  // Warp images
  m_OutputImages.clear();
  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
  {
    typedef itk::WarpImageFilter<
      ImageType, ImageType, DeformationFieldType>
      WarperType;
    typename WarperType::Pointer warpf = WarperType::New();
    warpf->SetInput(m_MovingImages[ichan]);
    warpf->SetDeformationField(dispF);
    warpf->SetOutputDirection(m_DeformationField->GetDirection());
    warpf->SetOutputOrigin(m_DeformationField->GetOrigin());
    warpf->SetOutputSpacing(m_DeformationField->GetSpacing());
    warpf->Update();
    m_OutputImages.push_back(warpf->GetOutput());
  }
}

#endif