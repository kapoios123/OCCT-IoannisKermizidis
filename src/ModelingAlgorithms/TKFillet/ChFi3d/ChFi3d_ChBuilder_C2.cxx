// Created on: 1996-07-26
// Created by: Stagiaire Xuan Trang PHAMPHU
// Copyright (c) 1996-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#include <ChFi3d_ChBuilder.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <ChFi3d_Builder_0.hxx>
#include <ChFiDS_CommonPoint.hxx>
#include <ChFiDS_FaceInterference.hxx>
#include <ChFiDS_Regul.hxx>
#include <ChFiDS_Spine.hxx>
#include <ChFiDS_State.hxx>
#include <ChFiDS_Stripe.hxx>
#include <ChFiDS_SurfData.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomFill_Boundary.hxx>
#include <GeomFill_ConstrainedFilling.hxx>
#include <NCollection_HSequence.hxx>
#include <Standard_ErrorHandler.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs.hxx>
#include <TopExp_Explorer.hxx>
#include <TopOpeBRepDS_Curve.hxx>
#include <TopOpeBRepDS_DataStructure.hxx>
#include <TopOpeBRepDS_HDataStructure.hxx>
#include <TopOpeBRepDS_SurfaceCurveInterference.hxx>
#include <TopoDS.hxx>

//=================================================================================================

void ChFi3d_ChBuilder::PerformTwoCorner(const int Index)
{
  if (!PerformTwoCornerSameDiff(Index) && !PerformTwoCornerbyInter(Index))
  {
    PerformMoreThreeCorner(Index, 2);
    done = true;
  }
}

//=================================================================================================

bool ChFi3d_ChBuilder::PerformTwoCornerSameDiff(const int Index)
{
  const TopoDS_Vertex&        Vtx  = myVDataMap.FindKey(Index);
  TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();

  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It(myVDataMap(Index));
  occ::handle<ChFiDS_Stripe>                             st1 = It.Value();
  It.Next();
  occ::handle<ChFiDS_Stripe> st2 = It.Value();
  if (st1 == st2)
  {
    return false;
  }
  int        Sens1 = 0, Sens2 = 0;
  const int  Isd1     = ChFi3d_IndexOfSurfData(Vtx, st1, Sens1);
  const int  Isd2     = ChFi3d_IndexOfSurfData(Vtx, st2, Sens2);
  const bool isfirst1 = (Sens1 == 1);
  const bool isfirst2 = (Sens2 == 1);

  const ChFiDS_State Stat1 = st1->Spine()->Status(isfirst1);
  const ChFiDS_State Stat2 = st2->Spine()->Status(isfirst2);
  bool               samefirst;
  if (Stat1 == ChFiDS_OnSame && Stat2 == ChFiDS_OnDiff)
  {
    samefirst = true;
  }
  else if (Stat1 == ChFiDS_OnDiff && Stat2 == ChFiDS_OnSame)
  {
    samefirst = false;
  }
  else
  {
    return false;
  }

  bool        Okvisavis = false, SameSide = false;
  int         IFaCo1 = 0, IFaCo2 = 0;
  double      UIntPC1 = 0., UIntPC2 = 0.;
  TopoDS_Face FaCo;
  const bool  OkinterCC = ChFi3d_IsInFront(DStr,
                                          st1,
                                          st2,
                                          Isd1,
                                          Isd2,
                                          Sens1,
                                          Sens2,
                                          UIntPC1,
                                          UIntPC2,
                                          FaCo,
                                          SameSide,
                                          IFaCo1,
                                          IFaCo2,
                                          Okvisavis,
                                          Vtx,
                                          true);
  if (!Okvisavis || !OkinterCC || SameSide)
  {
    return false;
  }

  const occ::handle<ChFiDS_Stripe>& stsam = samefirst ? st1 : st2;
  const occ::handle<ChFiDS_Stripe>& stdif = samefirst ? st2 : st1;
  occ::handle<ChFiDS_SurfData>      sdsam =
    stsam->SetOfSurfData()->Value(samefirst ? Isd1 : Isd2);
  occ::handle<ChFiDS_SurfData> sddif = stdif->SetOfSurfData()->Value(samefirst ? Isd2 : Isd1);
  const double                 uintpcsam  = samefirst ? UIntPC1 : UIntPC2;
  const double                 uintpcdif  = samefirst ? UIntPC2 : UIntPC1;
  const int                    ifacosam   = samefirst ? IFaCo1 : IFaCo2;
  const int                    ifacodif   = samefirst ? IFaCo2 : IFaCo1;
  const int                    ifaopsam   = 3 - ifacosam;
  const int                    ifaopdif   = 3 - ifacodif;
  const bool                   isfirstsam = samefirst ? isfirst1 : isfirst2;
  const bool                   isfirstdif = samefirst ? isfirst2 : isfirst1;

  ChFiDS_CommonPoint& cpopdif = sddif->ChangeVertex(isfirstdif, ifaopdif);
  if (!cpopdif.IsOnArc())
  {
    return false;
  }
  const TopoDS_Edge& Arcopdif = cpopdif.Arc();
  const TopoDS_Face& Fopsam   = TopoDS::Face(DStr.Shape(sdsam->Index(ifaopsam)));
  TopExp_Explorer    ex(Fopsam, TopAbs_EDGE);
  for (; ex.More() && !ex.Current().IsSame(Arcopdif); ex.Next())
  {
  }
  if (!ex.More())
  {
    return false;
  }

  occ::handle<Geom_Surface>        surfsam  = DStr.Surface(sdsam->Surf()).Surface();
  occ::handle<Geom_Surface>        surfdif  = DStr.Surface(sddif->Surf()).Surface();
  occ::handle<GeomAdaptor_Surface> Hsurfsam = new GeomAdaptor_Surface(surfsam);
  occ::handle<GeomAdaptor_Surface> Hsurfdif = new GeomAdaptor_Surface(surfdif);
  occ::handle<BRepAdaptor_Surface> HBRFopsam = new BRepAdaptor_Surface(Fopsam, false);
  occ::handle<Geom2d_Curve>        pcsurfsam, pcsurfdif, pcFopsam;
  occ::handle<Geom_Surface>        Surfcoin;
  try
  {
    OCC_CATCH_SIGNALS
    const gp_Pnt2d ppopsam = sdsam->Interference(ifaopsam).PCurveOnSurf()->Value(uintpcsam);
    const gp_Pnt2d ppcosam = sdsam->Interference(ifacosam).PCurveOnSurf()->Value(uintpcsam);
    occ::handle<GeomFill_Boundary> Bsam =
      ChFi3d_mkbound(Hsurfsam, pcsurfsam, ppopsam, ppcosam, tolapp3d, 2.e-4, true);
    const double   upcopdif = sddif->Interference(ifaopdif).Parameter(isfirstdif);
    const gp_Pnt2d ppopdif  = sddif->Interference(ifaopdif).PCurveOnSurf()->Value(upcopdif);
    const gp_Pnt2d ppcodif  = sddif->Interference(ifacodif).PCurveOnSurf()->Value(uintpcdif);
    occ::handle<GeomFill_Boundary> Bdif =
      ChFi3d_mkbound(Hsurfdif, pcsurfdif, ppcodif, ppopdif, tolapp3d, 2.e-4, true);

    gp_Pnt2d ppfacsam, ppfacdif;
    sdsam->Interference(ifaopsam).PCurveOnFace()->D0(uintpcsam, ppfacsam);
    BRepAdaptor_Curve2d PCArcFac(Arcopdif, Fopsam);
    PCArcFac.D0(cpopdif.ParameterOnArc(), ppfacdif);
    if (HBRFopsam->IsUClosed())
    {
      const double Uperiod = HBRFopsam->LastUParameter() - HBRFopsam->FirstUParameter();
      if (std::abs(ppfacsam.X() - ppfacdif.X()) > Uperiod / 2)
      {
        ppfacdif.SetX(ppfacdif.X() + (ppfacdif.X() < ppfacsam.X() ? Uperiod : -Uperiod));
      }
    }
    occ::handle<GeomFill_Boundary> Bfac =
      ChFi3d_mkbound(HBRFopsam, pcFopsam, ppfacsam, ppfacdif, tolapp3d, 2.e-4, true);

    GeomFill_ConstrainedFilling fil(8, 20);
    fil.Init(Bsam, Bdif, Bfac, true);
    Surfcoin = fil.Surface();
  }
  catch (Standard_Failure const&)
  {
    return false;
  }
  if (Surfcoin.IsNull())
  {
    return false;
  }

  occ::handle<ChFiDS_Stripe> corner = new ChFiDS_Stripe();
  occ::handle<NCollection_HSequence<occ::handle<ChFiDS_SurfData>>>& cornerset =
    corner->ChangeSetOfSurfData();
  cornerset                         = new NCollection_HSequence<occ::handle<ChFiDS_SurfData>>();
  occ::handle<ChFiDS_SurfData> coin = new ChFiDS_SurfData();
  cornerset->Append(coin);

  occ::handle<Geom2d_Curve> pcnul;
  if (!CompleteData(coin,
                    Surfcoin,
                    Hsurfsam,
                    pcsurfsam,
                    HBRFopsam,
                    pcnul,
                    sdsam->Orientation(),
                    true,
                    false,
                    false,
                    false,
                    false))
  {
    throw Standard_Failure("TwoCorner : chamfer corner filling failed");
  }

  double                    P1deb, P2deb, P1fin, P2fin, tolreached;
  const ChFiDS_CommonPoint& Pf1 = coin->VertexFirstOnS1();
  ChFiDS_CommonPoint&       Pf2 = coin->ChangeVertexFirstOnS2();
  const ChFiDS_CommonPoint& Pl1 = coin->VertexLastOnS1();
  ChFiDS_CommonPoint&       Pl2 = coin->ChangeVertexLastOnS2();
  Pf2 = Pl2    = cpopdif;
  const int If1 = ChFi3d_IndexPointInDS(Pf1, DStr);
  const int If2 = ChFi3d_IndexPointInDS(Pf2, DStr);
  const int Il1 = ChFi3d_IndexPointInDS(Pl1, DStr);
  const int Il2 = If2;

  gp_Pnt2d pp1 =
    coin->InterferenceOnS1().PCurveOnSurf()->Value(coin->InterferenceOnS1().FirstParameter());
  gp_Pnt2d pp2 =
    coin->InterferenceOnS2().PCurveOnSurf()->Value(coin->InterferenceOnS2().FirstParameter());
  occ::handle<Geom_Curve> C3d;
  ChFi3d_ComputeArete(Pf1,
                      pp1,
                      Pf2,
                      pp2,
                      DStr.Surface(coin->Surf()).Surface(),
                      C3d,
                      corner->ChangeFirstPCurve(),
                      P1deb,
                      P2deb,
                      tolapp3d,
                      tol2d,
                      tolreached,
                      0);
  double                         tolr1;
  occ::handle<GeomAdaptor_Curve> HC3d = new GeomAdaptor_Curve(C3d);
  ChFi3d_SameParameter(HC3d, pcFopsam, HBRFopsam, tolapp3d, tolr1);
  tolreached = std::max(tolreached, tolr1);
  const int          Icf = DStr.AddCurve(TopOpeBRepDS_Curve(C3d, tolreached));
  TopAbs_Orientation OpcFopsam = sdsam->Interference(ifaopsam).Transition();
  const int          IFopsam   = sdsam->Index(ifaopsam);
  if (isfirstsam)
  {
    OpcFopsam = TopAbs::Reverse(OpcFopsam);
  }
  DStr.ChangeShapeInterferences(IFopsam).Append(
    ChFi3d_FilCurveInDS(Icf, IFopsam, pcFopsam, OpcFopsam));

  ChFiDS_Regul regdeb, regfin;
  regdeb.SetCurve(Icf);
  regdeb.SetS1(coin->Surf(), false);
  regdeb.SetS2(IFopsam, true);
  myRegul.Append(regdeb);
  corner->ChangeFirstCurve(Icf);
  corner->ChangeFirstParameters(P1deb, P2deb);
  corner->ChangeIndexFirstPointOnS1(If1);
  corner->ChangeIndexFirstPointOnS2(If2);

  pp1 = coin->InterferenceOnS1().PCurveOnSurf()->Value(coin->InterferenceOnS1().LastParameter());
  pp2 = coin->InterferenceOnS2().PCurveOnSurf()->Value(coin->InterferenceOnS2().LastParameter());
  ChFi3d_ComputeArete(Pl1,
                      pp1,
                      Pl2,
                      pp2,
                      DStr.Surface(coin->Surf()).Surface(),
                      C3d,
                      corner->ChangeLastPCurve(),
                      P1fin,
                      P2fin,
                      tolapp3d,
                      tol2d,
                      tolreached,
                      0);
  double tolr2;
  HC3d->Load(C3d);
  ChFi3d_SameParameter(HC3d, pcsurfdif, Hsurfdif, tolapp3d, tolr2);
  tolreached    = std::max(tolreached, tolr2);
  const int Icl = DStr.AddCurve(TopOpeBRepDS_Curve(C3d, tolreached));
  regfin.SetCurve(Icl);
  regfin.SetS1(coin->Surf(), false);
  regfin.SetS2(sddif->Surf(), false);
  myRegul.Append(regfin);
  corner->ChangeLastCurve(Icl);
  corner->ChangeLastParameters(P1fin, P2fin);
  corner->ChangeIndexLastPointOnS1(Il1);
  corner->ChangeIndexLastPointOnS2(Il2);

  coin->ChangeIndexOfS1(-sdsam->Surf());
  coin->ChangeIndexOfS2(0);
  corner->SetSolidIndex(stsam->SolidIndex());

  const ChFiDS_FaceInterference& intcoin1 = coin->InterferenceOnS1();
  stsam->SetCurve(intcoin1.LineIndex(), isfirstsam);
  stsam->InDS(isfirstsam);
  stsam->ChangePCurve(isfirstsam) = intcoin1.PCurveOnFace();
  stsam->SetIndexPoint(If1, isfirstsam, ifaopsam);
  stsam->SetIndexPoint(Il1, isfirstsam, ifacosam);
  stsam->SetParameters(isfirstsam, intcoin1.FirstParameter(), intcoin1.LastParameter());
  sdsam->ChangeVertex(isfirstsam, ifaopsam) = Pf1;
  sdsam->ChangeVertex(isfirstsam, ifacosam) = Pl1;
  sdsam->ChangeInterferenceOnS1().SetParameter(uintpcsam, isfirstsam);
  sdsam->ChangeInterferenceOnS2().SetParameter(uintpcsam, isfirstsam);
  if (ifaopsam == 2)
  {
    stsam->SetOrientation(TopAbs_REVERSED, isfirstsam);
  }

  stdif->SetCurve(Icl, isfirstdif);
  stdif->ChangePCurve(isfirstdif) = pcsurfdif;
  stdif->SetIndexPoint(Il2, isfirstdif, ifaopdif);
  stdif->SetIndexPoint(Il1, isfirstdif, ifacodif);
  stdif->SetParameters(isfirstdif, P1fin, P2fin);
  sddif->ChangeVertex(isfirstdif, ifaopdif) = Pl2;
  sddif->ChangeVertex(isfirstdif, ifacodif) = Pl1;
  sddif->ChangeInterference(ifacodif).SetParameter(uintpcdif, isfirstdif);
  if (ifaopdif == 1)
  {
    stdif->SetOrientation(TopAbs_REVERSED, isfirstdif);
  }

  if (!myEVIMap.IsBound(Vtx))
  {
    myEVIMap.Bind(Vtx, NCollection_List<int>());
  }
  myEVIMap.ChangeFind(Vtx).Append(coin->Surf());
  myListStripe.Append(corner);
  done = true;
  return true;
}
