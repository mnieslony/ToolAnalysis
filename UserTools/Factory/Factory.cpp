#include "Factory.h"

Tool* Factory(std::string tool) {
Tool* ret=0;

// if (tool=="Type") tool=new Type;
if (tool=="DummyTool") ret=new DummyTool;
if (tool=="ExampleGenerateData") ret=new ExampleGenerateData;
if (tool=="ExampleSaveStore") ret=new ExampleSaveStore;
if (tool=="ExampleSaveRoot") ret=new ExampleSaveRoot;
if (tool=="ExampleloadStore") ret=new ExampleloadStore;
if (tool=="ExamplePrintData") ret=new ExamplePrintData;
if (tool=="ExampleLoadRoot") ret=new ExampleLoadRoot;
if (tool=="PythonScript") ret=new PythonScript;
if (tool=="LAPPDParseScope") ret=new LAPPDParseScope;
if (tool=="LAPPDParseACC") ret=new LAPPDParseACC;
if (tool=="LAPPDFindPeak") ret=new LAPPDFindPeak;
if (tool=="SaveANNIEEvent") ret=new SaveANNIEEvent;
if (tool=="LAPPDSim") ret=new LAPPDSim;
if (tool=="LoadWCSim") ret=new LoadWCSim;
if (tool=="FindMrdTracks") ret=new FindMrdTracks;
if (tool=="PrintANNIEEvent") ret=new PrintANNIEEvent;
if (tool=="GenerateHits") ret=new GenerateHits;
if (tool=="LAPPDcfd") ret=new LAPPDcfd;
if (tool=="LAPPDBaselineSubtract") ret=new LAPPDBaselineSubtract;
if (tool=="NeutronStudyReadSandbox") ret=new NeutronStudyReadSandbox;
if (tool=="NeutronStudyPMCS") ret=new NeutronStudyPMCS;
if (tool=="NeutronStudyWriteTree") ret=new NeutronStudyWriteTree;
if (tool=="BeamTimeAna") ret=new BeamTimeAna;
if (tool=="BeamTimeTreeMaker") ret=new BeamTimeTreeMaker;
if (tool=="BeamTimeTreeReader") ret=new BeamTimeTreeReader;
if (tool=="RawLoader") ret=new RawLoader;
if (tool=="LAPPDBaselineSubtract") ret=new LAPPDBaselineSubtract;
if (tool=="LAPPDSaveROOT") ret=new LAPPDSaveROOT;
if (tool=="LAPPDFilter") ret=new LAPPDFilter;
if (tool=="LAPPDIntegratePulse") ret=new LAPPDIntegratePulse;
if (tool=="ADCCalibrator") ret=new ADCCalibrator;
if (tool=="ADCHitFinder") ret=new ADCHitFinder;
if (tool=="BeamChecker") ret=new BeamChecker;
if (tool=="BeamFetcher") ret=new BeamFetcher;
if (tool=="FindTrackLengthInWater") ret=new FindTrackLengthInWater;
if (tool=="LoadANNIEEvent") ret=new LoadANNIEEvent;
if (tool=="PhaseITreeMaker") ret=new PhaseITreeMaker;
if (tool=="MrdPaddlePlot") ret=new MrdPaddlePlot;
if (tool=="LoadWCSimLAPPD") ret=new LoadWCSimLAPPD;
if (tool=="WCSimDemo") ret=new WCSimDemo;
if (tool=="DigitBuilder") ret=new DigitBuilder;
if (tool=="VtxSeedGenerator") ret=new VtxSeedGenerator;
if (tool=="VtxPointPositionFinder") ret=new VtxPointPositionFinder;
if (tool=="LAPPDlasertestHitFinder") ret=new LAPPDlasertestHitFinder;
if (tool=="RawLoadToRoot") ret=new RawLoadToRoot;
if (tool=="MRDPulseFinder") ret=new MRDPulseFinder;
if (tool=="LAPPDAnalysis") ret=new LAPPDAnalysis;
if (tool=="ExampleOverTool") ret=new ExampleOverTool;
if (tool=="PhaseIITreeMaker") ret=new PhaseIITreeMaker;
if (tool=="VertexGeometryCheck") ret=new VertexGeometryCheck;
if (tool=="LikelihoodFitterCheck") ret=new LikelihoodFitterCheck;
if (tool=="EventSelector") ret=new EventSelector;
if (tool=="SaveRecoEvent") ret=new SaveRecoEvent;
if (tool=="VtxExtendedVertexFinder") ret=new VtxExtendedVertexFinder;
if (tool=="VtxPointDirectionFinder") ret=new VtxPointDirectionFinder;
if (tool=="VtxPointVertexFinder") ret=new VtxPointVertexFinder;
if (tool=="LoadCCData") ret=new LoadCCData;
if (tool=="WaveformNNLS") ret=new WaveformNNLS;
if (tool=="HitCleaner") ret=new HitCleaner;
if (tool=="HitResiduals") ret=new HitResiduals;
if (tool=="MonitorReceive") ret=new MonitorReceive;
if (tool=="MonitorSimReceive") ret=new MonitorSimReceive;
if (tool=="DigitBuilderDoE") ret=new DigitBuilderDoE;
if (tool=="EventSelectorDoE") ret=new EventSelectorDoE;
if (tool=="MonitorMRDTime") ret=new MonitorMRDTime;
if (tool=="MonitorMRDLive") ret=new MonitorMRDLive;
if (tool=="PulseSimulation") ret=new PulseSimulation;
if (tool=="PlotLAPPDTimesFromStore") ret=new PlotLAPPDTimesFromStore;
if (tool=="CheckDetectorCounts") ret=new CheckDetectorCounts;
if (tool=="MrdDistributions") ret=new MrdDistributions;
if (tool=="MCParticleProperties") ret=new MCParticleProperties;
if (tool=="DigitBuilderROOT") ret=new DigitBuilderROOT;
if (tool=="MrdEfficiency") ret=new MrdEfficiency;
if (tool=="EventDisplay") ret=new EventDisplay;
if (tool=="TankCalibrationDiffuser") ret=new TankCalibrationDiffuser;
if (tool=="TotalLightMap") ret=new TotalLightMap;
if (tool=="MrdDiscriminatorScan") ret=new MrdDiscriminatorScan;
if (tool=="MCRecoEventLoader") ret=new MCRecoEventLoader;
if (tool=="MonitorMRDEventDisplay") ret=new MonitorMRDEventDisplay;
if (tool=="LoadGeometry") ret=new LoadGeometry;
if (tool=="LoadRATPAC") ret=new LoadRATPAC;
if (tool=="TimeClustering") ret=new TimeClustering;
if (tool=="GracefulStop") ret=new GracefulStop;
if (tool=="PhaseIIADCHitFinder") ret=new PhaseIIADCHitFinder;
if (tool=="TrackCombiner") ret=new TrackCombiner;
if (tool=="SimulatedWaveformDemo") ret=new SimulatedWaveformDemo;
if (tool=="CNNImage") ret=new CNNImage;
if (tool=="MonitorTankTime") ret=new MonitorTankTime;
if (tool=="PhaseIIADCCalibrator") ret=new PhaseIIADCCalibrator;
if (tool=="MCHitToHitComparer") ret=new MCHitToHitComparer;
if (tool=="MCPropertiesToTree") ret=new MCPropertiesToTree;
if (tool=="CalcClassificationVars") ret=new CalcClassificationVars;
if (tool=="StoreClassificationVars") ret=new StoreClassificationVars;
if (tool=="LoadGenieEvent") ret=new LoadGenieEvent;
if (tool=="PrintGenieEvent") ret=new PrintGenieEvent;
if (tool=="PlotWaveforms") ret=new PlotWaveforms;
if (tool=="PMTDataDecoder") ret=new PMTDataDecoder;
if (tool=="ANNIEEventBuilder") ret=new ANNIEEventBuilder;
if (tool=="MRDDataDecoder") ret=new MRDDataDecoder;
if (tool=="PrintADCData") ret=new PrintADCData;
if (tool=="ClusterFinder") ret=new ClusterFinder;
if (tool=="PrintRecoEvent") ret=new PrintRecoEvent;
if (tool=="RunValidation") ret=new RunValidation;
if (tool=="AmBeRunStatistics") ret=new AmBeRunStatistics;
if (tool=="SimpleTankEnergyCalibrator") ret=new SimpleTankEnergyCalibrator;
if (tool=="BeamClusterPlots") ret=new BeamClusterPlots;
if (tool=="MrdPaddleEfficiencyPreparer") ret=new MrdPaddleEfficiencyPreparer;
if (tool=="MrdPaddleEfficiencyCalc") ret=new MrdPaddleEfficiencyCalc;
if (tool=="FMVEfficiency") ret=new FMVEfficiency;
if (tool=="LoadRawData") ret=new LoadRawData;
if (tool=="TriggerDataDecoder") ret=new TriggerDataDecoder;
if (tool=="ClusterClassifiers") ret=new ClusterClassifiers;
if (tool=="MRDLoopbackAnalysis") ret=new MRDLoopbackAnalysis;
if (tool=="VetoEfficiency") ret=new VetoEfficiency;
if (tool=="MonitorTrigger") ret=new MonitorTrigger;

if (tool=="EventClassification") ret=new EventClassification;

if (tool=="DataSummary") ret=new DataSummary;
if (tool=="MaxPEPlots") ret=new MaxPEPlots;
if (tool=="PlotDecodedTimestamps") ret=new PlotDecodedTimestamps;
if (tool=="StoreDecodedTimestamps") ret=new StoreDecodedTimestamps;
return ret;
}
