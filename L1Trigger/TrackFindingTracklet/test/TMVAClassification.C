void TMVAClassification( TString myMethodList = "" )
{
   // This loads the library
   TMVA::Tools::Instance();

   // This loads the input file for machine learning
   TFile * input = TFile::Open("output_TTbar_PU200_WithTruncation_12dec2018.root", "read");
   
   // Register the training and test trees
   // signal     -> electron
   // background -> muon
   TTree *background = (TTree*)input->Get("t_fake");
   TTree *signal     = (TTree*)input->Get("t_real");

   // Create a output file where TMVA results will be stored 
   TString outputFileName = "tmva_output1.root";
   TFile* output = TFile::Open( outputFileName, "recreate" );

   // Create the factory object where methods for machine learning are contained
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification",
					       output,
					       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // Create the dataloader where train/test objects is handled
   TMVA::DataLoader *dataloader=new TMVA::DataLoader("MLP");

   // Add varaible for machine learning in Float ('F') type
   dataloader->AddVariable( "pt := pt", 'F' );
   dataloader->AddVariable( "eta := eta", 'F' );
   dataloader->AddVariable( "nstub := nstub", 'I' );
   dataloader->AddVariable( "chi2 := chi2", 'F' );

   // Set event weights: Just make it 1 
   Float_t signalWeight     = 1.0;
   Float_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signal    ,     signalWeight );
   dataloader->AddBackgroundTree( background, backgroundWeight );

   dataloader->PrepareTrainingAndTestTree( "","",
					   "nTrain_Signal=5000:nTrain_Background=5000:nTest_Signal=5000:nTest_Background=5000:SplitMode=Random:NormMode=NumEvents:!V" );

   // ### Book MVA methods
   // Every parameter for machine learning is roughly optimized already so you don't have to edit
//   factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP4", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N-1:TestRate=5:!UseRegulator" );
//   factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP8", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+3:TestRate=5:!UseRegulator" );
//   factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP12", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+7:TestRate=5:!UseRegulator" );
////   factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT5",
////            "!H:!V:NTrees=5:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
////   factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT850",
////                        "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
//   factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT2500",
//                        "!H:!V:NTrees=2500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
    factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM0.01", "Gamma=0.01:C=1:Tol=0.01:MaxIter=1000" );
    factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM0.03", "Gamma=0.03:C=1:Tol=0.01:MaxIter=1000" );//the best
   factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM1", "Gamma=1:C=1:Tol=0.01:MaxIter=1000" );
////    // deep network
////   factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN1", "Layout=TANH|3,LINEAR");
////   factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN1", "Layout=TANH|10,LINEAR");
////   factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN1", "Layout=TANH|10,TANH|10,TANH|10,LINEAR");

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();
    
   // Save the output
   output->Close();

   delete factory;
   delete dataloader;

   // Pop-up GUI of TMVA result
   TMVA::TMVAGui( outputFileName );
}
