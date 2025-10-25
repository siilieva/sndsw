void UserDecayConfig() {
   cout << "Loading User Decay Config from macro"<< endl;  
   TDatabasePDG *db= TDatabasePDG::Instance();
   TParticlePDG *p=0;

   Int_t mode[6][3];                  
   Float_t bratio[6];
   Int_t DplusPDG;
   p= db->GetParticle("D+"); //411
   if(p) DplusPDG=p->PdgCode();// why not just set the correct pdg code?
   for (Int_t kz = 0; kz < 6; kz++) {
      bratio[kz] = 0.;
      mode[kz][0] = 0;
      mode[kz][1] = 0;
      mode[kz][2] = 0;
    //  cout << mode[kz][0] << " " << 	mode[kz][1] << " " << mode[kz][2] << endl;
   }
   bratio[0] = 100.;// in percent
   mode[0][0] =-13  ;
   mode[0][1] =14  ;
   
 /*  bratio[1] = 50.;
   mode[1][0] =2212  ;
   mode[1][1] =DplusPDG  ;
    
  */
   gMC->SetDecayMode(DplusPDG,bratio,mode);
    

    
}
