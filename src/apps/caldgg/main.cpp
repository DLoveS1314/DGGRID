#include <iostream>

using namespace std;

#include <dglib/DgIDGGS4H.h>
#include <dglib/DgIDGGS4D.h>
#include <dglib/DgIDGG.h>
#include "dglib/DgBoundedIDGG.h"
#include "dglib/DgBoundedIDGGS.h"
#include "calDgg.h"
#include "testfun.h"

#include <omp.h>
#include <iostream>

#include <dglib/DgApSeq.h>
vector<string> inFilepaths ={
    "/home/dls/data/openmmlab/mmclassification/tools/pandas_save/hascode/FULLER4D_6",
    "/home/dls/data/openmmlab/mmclassification/tools/pandas_save/hascode/FULLER4H_6",
    "/home/dls/data/openmmlab/mmclassification/tools/pandas_save/hascode/FULLER4T_6",

    "/home/dls/data/openmmlab/mmclassification/tools/pandas_save/hascode/ISEA4D_6",
    "/home/dls/data/openmmlab/mmclassification/tools/pandas_save/hascode/ISEA4H_6",
    "/home/dls/data/openmmlab/mmclassification/tools/pandas_save/hascode/ISEA4T_6"
};
vector<string> outFileNames ={
    "/home/dls/data/openmmlab/DGGRID/src/apps/caldgg/c/FULLER4D_6_caldgg.csv",
    "/home/dls/data/openmmlab/DGGRID/src/apps/caldgg/c/FULLER4H_6_caldgg.csv",
    "/home/dls/data/openmmlab/DGGRID/src/apps/caldgg/c/FULLER4T_6_caldgg.csv",

    "/home/dls/data/openmmlab/DGGRID/src/apps/caldgg/c/ISEA4D_6_caldgg.csv",
    "/home/dls/data/openmmlab/DGGRID/src/apps/caldgg/c/ISEA4H_6_caldgg.csv",
    "/home/dls/data/openmmlab/DGGRID/src/apps/caldgg/c/ISEA4T_6_caldgg.csv"
};
void run(int aperture,DgGridTopology& gridTopo ,string  name ,string projType ,string inFilepath,string outFileName, int ptsPerEdgeDensify)
{
    
    // DgGridMetric gridMetricD = dgg::topo::D4;
    // string dggtype ="ISEA4T"
    bool isMixed43=false;
    int numAp4 = 0;
    bool isApSeq =false;
    DgApSeq apSeq;
    bool isSuperfund =false;
    DgRFNetwork net0;
    const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, "GS0"));

    DgGeoCoord vert0(0.0L, 90.0L, false); // args: lon, lat, isRadians
    long double azimuth = 0.0L;

    string suffix ="csv";
    string inaddtype= "SEQNUM";
    string outAddType = "Q2DI";
    char outputDelimiter =',';
    int res =6;
    DgGridMetric gridMetric = dgg::topo::D4; //随便传只有菱形需要用
    
    
    const DgIDGGSBase *dggs  = DgIDGGSBase::makeRF(net0, geoRF,  vert0,
             azimuth,  aperture, res+3,  gridTopo,
             gridMetric, name,  projType,  isMixed43,  numAp4,
              isSuperfund,  isApSeq,  apSeq);
    const DgIDGGBase& dgg = dggs ->idggBase(res);

    calDgg cdgg = calDgg(dgg);
    cdgg.result2file1(outFileName, outputDelimiter,inFilepath, inaddtype,outAddType,suffix,ptsPerEdgeDensify)  ;
}
void calresult()
{

    int i = 0;
    string name ="FULLER4D";
    string projType = "FULLER";
    int aperture =4;
    DgGridTopology gridTopo =DgGridTopology::Diamond;
    std::string inFilepath = inFilepaths[i]; // 替换为实际的文件夹路径
    string outFileName = outFileNames[i];
    int ptsPerEdgeDensify =0;

    cout<<name<<","<<projType<<","<<inFilepath<<","<<outFileName<<","<<ptsPerEdgeDensify<<endl;

    // run( aperture,  gridTopo ,   name ,  projType , inFilepath,  outFileName,ptsPerEdgeDensify);

    cout<<"======"<<name<<" complete======"<<endl;

    i = i+1;//1
    name ="FULLER4H";
    projType = "FULLER";
    aperture =4;
    gridTopo =DgGridTopology::Hexagon;
    inFilepath = inFilepaths[i]; // 替换为实际的文件夹路径
    outFileName = outFileNames[i];
     ptsPerEdgeDensify =0;

    cout<<name<<","<<projType<<","<<inFilepath<<","<<outFileName<<","<<ptsPerEdgeDensify<<endl;


    // run( aperture,  gridTopo ,   name ,  projType , inFilepath,  outFileName,ptsPerEdgeDensify);
 
    cout<<"======"<<name<<" complete======"<<endl;


    i = i+1;//2
    name ="FULLER4T";
    projType = "FULLER";
    aperture =4;
    gridTopo =DgGridTopology::Triangle;
    inFilepath = inFilepaths[i]; // 替换为实际的文件夹路径
    outFileName = outFileNames[i];
      ptsPerEdgeDensify =0;

    cout<<name<<","<<projType<<","<<inFilepath<<","<<outFileName<<","<<ptsPerEdgeDensify<<endl;


    // run( aperture,  gridTopo ,   name ,  projType , inFilepath,  outFileName,ptsPerEdgeDensify);
    cout<<"======"<<name<<" complete======"<<endl;
 


    i = i+1;//3
    name ="ISEA4D";
    projType = "ISEA";
    aperture =4;
    gridTopo =DgGridTopology::Diamond;
    inFilepath = inFilepaths[i]; // 替换为实际的文件夹路径
    outFileName = outFileNames[i];
     ptsPerEdgeDensify =50;
    
    cout<<name<<","<<projType<<","<<inFilepath<<","<<outFileName<<","<<ptsPerEdgeDensify<<endl;


    // run( aperture,  gridTopo ,   name ,  projType , inFilepath,  outFileName,ptsPerEdgeDensify);
    cout<<"======"<<name<<" complete======"<<endl;
 


    i = i+1;//4
    name ="ISEA4H";
    projType = "ISEA";
    aperture =4;
    gridTopo =DgGridTopology::Hexagon;
    inFilepath = inFilepaths[i]; // 替换为实际的文件夹路径
    outFileName = outFileNames[i];
    ptsPerEdgeDensify =50;

    cout<<name<<","<<projType<<","<<inFilepath<<","<<outFileName<<","<<ptsPerEdgeDensify<<endl;


    // run( aperture,  gridTopo ,   name ,  projType , inFilepath,  outFileName,ptsPerEdgeDensify);
    cout<<"======"<<name<<" complete======"<<endl;
 


    i = i+1;//5
    name ="ISEA4T";
    projType = "ISEA";
    aperture =4;
    gridTopo =DgGridTopology::Triangle;
    inFilepath = inFilepaths[i]; // 替换为实际的文件夹路径
    outFileName = outFileNames[i];
      ptsPerEdgeDensify =50;

    cout<<name<<","<<projType<<","<<inFilepath<<","<<outFileName<<","<<ptsPerEdgeDensify<<endl;


    run( aperture,  gridTopo ,   name ,  projType , inFilepath,  outFileName,ptsPerEdgeDensify);
    cout<<"======"<<name<<" complete======"<<endl;
 

  


}
void testomp()
{
    int num_threads = 4;

    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();

        #pragma omp for
        for (int i = 0; i < 100; i++) {
            std::cout << "Thread " << thread_id << " of " << num_threads << " is processing index " << i << std::endl;
        }
    }


}
 

int main (int, char**)
{
//    testIDGGbounds();
    // testcalArea();
    // testcalange();
    // testreadadd();
    // testfileio();
    // testtrans2QDI();
    // testiteratefile();
    // testresult2file();
    // testomp();
    calresult();
    // printf("end");
    // DgRFNetwork net0;
    // const DgGeoSphRF& geoRF = *(DgGeoSphRF::makeRF(net0, "GS0"));
    // DgGeoCoord vert0(0.0L, 90.0L, false); // args: lon, lat, isRadians

    // long double azimuth = 0.0L;
    // const DgIDGGS4D * idggsDPtr = DgIDGGS4D::makeRF(net0, geoRF, vert0, azimuth, 10,"ISEA4D_D4","ISEA",dgg::topo::D4);
    // const DgIDGGS4D& idggs = *idggsDPtr;
    // const DgIDGGBase& dgg = idggs.idgg(6);


    // dgcerr<<dgg.gridTopo()<< " setAddNeighbors not yet implemented"<<endl;

    // dgcerr<<dgg.gridTopo()<< " setAddNeighbors not yet implemented"<<endl;
}


