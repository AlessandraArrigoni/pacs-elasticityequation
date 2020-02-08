#include "../include/Problem.h"

Problem::Problem (const GetPot& dataFile, Bulk* bulk1, Bulk* bulk2):
			M_Bulk1(bulk1),	M_Bulk2(bulk2),
			M_BC1(dataFile, "laplacian1/", "bulkData/"), // change the dataFile accordingly and set it
			M_BC2(dataFile, "laplacian2/", "bulkData/"), // change the dataFile accordingly and set it
			interfaceIdx1( dataFile(std::string("bulkData/laplacian1/interfaceIdx").data(),1) ),// change the dataFile accordingly and set it: devo capire come leggere il valore direttamente dal file!
			interfaceIdx2( dataFile(std::string("bulkData/laplacian2/interfaceIdx").data(),3)), // adesso le ho messe brutali così, poi andranno lette da file o altro.
			M_uFEM1( bulk1->getMesh(), dataFile, "femspaces/", "Sol", "bulkData/"),// change the dataFile accordingly and set it
			M_uFEM2( bulk2->getMesh(), dataFile, "femspaces/", "Sol", "bulkData/"),// change the dataFile accordingly and set it
			M_CoeffFEM1( bulk1->getMesh(), dataFile, "femspaces/", "Coeff", "bulkData/"), // change the dataFile accordingly and set it
			M_CoeffFEM2( bulk2->getMesh(), dataFile, "femspaces/", "Coeff", "bulkData/"), // change the dataFile accordingly and set it
			M_intMethod1(*(bulk1->getMesh()) ),
			M_intMethod2(*(bulk2->getMesh()) ),
		M_Sys()
{
	std::cout<<"l'indice dell'interfaccia su omega1 è "<<interfaceIdx1<<", quello su omega2 è "<<interfaceIdx2<<std::endl;

   M_nbDOF1 = M_uFEM1.nb_dof();
	 M_nbDOF2 = M_uFEM2.nb_dof();
	 M_nbTotDOF = M_nbDOF1 + M_nbDOF2;
	 std::cout<<"Il numero totale di dofs è "<<M_nbTotDOF<<std::endl;

	 Qdim = M_uFEM1.getFEM()->get_qdim(); // Dimensione dello spazio FEM (1 = scalare, 2 = vettoriale)

	 // Per ora faccio tutto sdoppiato, anche se poi probabilmente il metodo di integrazione sarà lo stesso sui 2 domini. Per ora ho solo copiato e incollato: in realtà la stringa non viene letta dal file di input ma prende il valore di default che comunque è lo stesso
  
  std::string intMethod(dataFile(std::string("bulkData/femspaces/integrationMethod").data(), "IM_TRIANGLE(19)" ) );
 std::cout<<"il metodo d'integrazione è:"<< intMethod<<std::endl;
   M_intMethod1.set_integration_method(bulk1->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod) );
	 M_intMethod2.set_integration_method(bulk2->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod) );

 
   M_Bulk1->getData()->setDiff(M_CoeffFEM1.getDOFpoints());  //valuta i coefficienti negli elementi della mesh
   M_Bulk2->getData()->setDiff(M_CoeffFEM2.getDOFpoints());  //valuta i coefficienti negli elementi della mesh


	 // Con questa funzione sto assegnando le facce della frontiera e dell'interfaccia a M_bulk1, M_bulk2 e alle sue regions quindi potrei cercare di recuperare il numero di dofs associati alla region che so essere l'interfaccia. Per adesso posso dire di sapere che è quella a destra per 1 e a sinistra per 2, poi magari posso scrivere una funzione che la determini in base alla valutazione dell'insieme di livello sui nodi.
   M_BC1.setBoundaries(M_Bulk1->getMesh());
	 M_BC2.setBoundaries(M_Bulk2->getMesh());

	 // Dato che per ipotesi le mesh sono conformi, il numero di dof sull'interfaccia è uguale per entrambi i domini, quindi per recuperarlo basta considerare solo uno dei due (a sinistra)

	 //Interface indices Omega1
 	size_type idxIFaceInDiriBD;
 	for (size_type j=0; j<M_BC1.getDiriBD().size(); j++){
 		if (M_BC1.getDiriBD()[j] == interfaceIdx1) {idxIFaceInDiriBD = j;}
 	}

	dal::bit_vector dal_dof_IFace1 = M_uFEM1.getFEM()->dof_on_region(M_BC1.getDiriBD()[idxIFaceInDiriBD]);
 	for(dal::bv_visitor i(dal_dof_IFace1); !i.finished(); ++i)
 	{		 M_rowsIFace1.push_back(i); 
			 M_rowsIFace.push_back(i);
	 }


	fromBitVectorToStdVector ( dal_dof_IFace1, dof_IFace1);
	std::cout<<"Gli indici sull'interfaccia di SINISTRA sono: "<<std::endl;
	for (size_type j=0; j<dof_IFace1.size(); j++){
	 std::cout<< dof_IFace1[j]<<"  ";
           }
	M_nbDOFIFace = dof_IFace1.size();

	std::cout<< "\nThe number of interface dofs is : "<<M_nbDOFIFace<<std::endl;
 	std::cout << "Interface Omega1      [OK]" << std::endl;

 	//Interface indices Omega2
 	size_type idxIFaceInNeumBD;
 	for (size_type j=0; j<M_BC2.getNeumBD().size(); j++){
 			if (M_BC2.getNeumBD()[j] == interfaceIdx2) {idxIFaceInNeumBD = j;}
 	}

 	dal::bit_vector dal_dof_IFace2 = M_uFEM2.getFEM()->dof_on_region(M_BC2.getNeumBD()[idxIFaceInNeumBD]);
 	for(dal::bv_visitor i(dal_dof_IFace2); !i.finished(); ++i)
 	{		M_rowsIFace2.push_back(i);	
			M_rowsIFace.push_back(i+M_nbDOF1); 
         }

	fromBitVectorToStdVector ( dal_dof_IFace2, dof_IFace2 );
	std::cout<<"Gli indici sull'interfaccia di DESTRA sono: "<<std::endl;
	for (size_type j=0; j<dof_IFace2.size(); j++){
	 std::cout<< dof_IFace2[j]<<"  ";
 }

 	std::cout << "Interface Omega2      [OK]" << std::endl;
}

void Problem::initialize()
{
         	M_uSol.reset(new scalarVector_Type (M_nbTotDOF));
         	gmm::clear(*M_uSol);

					M_uSol1.reset(new scalarVector_Type (M_nbDOF1));
         	gmm::clear(*M_uSol1);

					M_uSol2.reset(new scalarVector_Type (M_nbDOF2));
         	gmm::clear(*M_uSol2);
}

//restituisce un puntatore ai FEM del displacement.

FEM* Problem::getFEM(size_type const idx)
{
	if (idx == 1)
 		{return &M_uFEM1;}
	if (idx == 2)
		{return &M_uFEM2;}
	}

// dimensiona il sistema lineare per contenere il problema
void Problem::addToSys(LinearSystem* sys)
{
	M_Sys=sys;
        M_Sys->addToMatrix(M_nbTotDOF);

}

Bulk* Problem::getBulk(size_type const idx){
	if (idx == 1)
 		{return M_Bulk1;}
	if (idx == 2)
		{return M_Bulk2;}
}



// CONTO DEI DOF
// Il valore M_nbDOFIFace è quello relativo alle funzioni che vanno raddoppiate, quindi in questo caso M_nbTotDOF è dato dalla somma M_nbDOF1 + M_nbDOF2 per le trial functions visto che suppongo che siano 2 domini separati; poi per le test il numero di dofs sarà M_nbTotDOF - M_nbDOFIFace.
size_type Problem::getNDOF(std::string variable)
{
	if (variable=="u1")
	{		return M_nbDOF1;	}

	if (variable=="u2")
	{		return M_nbDOF2;	}

	if (variable=="v") // Meaning the test functions (smaller space)
	{		return M_nbTotDOF - M_nbDOFIFace;	}

	if (variable=="u") // Meaning the global solution (trial function)
	{		return M_nbTotDOF;	}

	if (variable=="all") // Equivalente a "u"
	{		return M_nbTotDOF; }
}

// assembla la matrice: decido di mettere prima tutti i dof relativi a Omega1 e poi quelli relativi a Omega2. Tra l'uno e l'altro (prima di imporre le condizioni di interfaccia) recupero le righe di Omega1 relative alle funzioni associate all'interfaccia perchè andranno spostate negli indici "globali" e sommate alle righe associate a Omega2.
void Problem::assembleMatrix(LinearSystem* sys)
{
	sparseMatrixPtr_Type A1, A2;

	A1.reset(new sparseMatrix_Type (M_nbDOF1, M_nbDOF1));
	gmm::clear(*A1);
	stiffness( A1, M_Bulk1, M_uFEM1, M_CoeffFEM1, M_intMethod1);
	
	A2.reset(new sparseMatrix_Type (M_nbDOF2, M_nbDOF2));
	gmm::clear(*A2);
	stiffness( A2, M_Bulk2, M_uFEM2, M_CoeffFEM2, M_intMethod2);





// dal blocco [A_1Gamma A_1GammaGamma 0 0] metto a zero i termini A_1GammaGamma e li sommo a A_2GammaGamma:
	for (size_type k=0; k<dof_IFace1.size(); k++)
	  {  
 	    size_type idx1=dof_IFace1[k];
             scalar_type value= (*A1)(idx1, idx1);
	   	size_type idx2=dof_IFace2[k];
                 // verifico che effettivamente i valori di Agammagamma1 sono uguali a quelli di Agammagamma2
		//std::cout<< "Valore di AGammaGamma1 è:" <<value<<std::endl;
	        //std::cout<< "Valore di AGammaGamma2 è:" <<(*A2)(idx2, idx2)<<std::endl;
	    (*A2)(idx2, idx2)+= value;
	   // (*A2)(idx2, idx2)/=2;
	    (*A1)(idx1, idx1)= 0;   // non serve mettere a zero questi termini perchè verranno comunque eliminati
              
	  }  
    

	M_Sys->addSubMatrix(A1, 0, 0);
	M_Sys->addSubMatrix(A2, M_nbDOF1, M_nbDOF1);
	
	

	// collego i termini legati all'interfaccia del dominio di sinistra a quelli dell'interfaccia del dominio di destra

	    for (size_type k=0; k<dof_IFace1.size(); k++)
		{
			sparseMatrixPtr_Type curRow1(new sparseMatrix_Type(1,M_nbTotDOF)); 
			sparseMatrixPtr_Type curCol2(new sparseMatrix_Type(M_nbTotDOF,1)); 
		         M_Sys->extractSubMatrix( dof_IFace1[k], 1, 0, M_nbTotDOF, curRow1 );     
			 M_Sys->extractSubMatrix( 0, M_nbTotDOF, dof_IFace1[k], 1, curCol2 );  
		         M_Sys->addSubMatrix(curRow1, dof_IFace2[k] + M_nbDOF1, 0);  
			 M_Sys->addSubMatrix(curCol2, 0, dof_IFace2[k]+M_nbDOF1);
			//M_Sys->setNullRow(dof_IFace1[k]);  // non serve mettere a zero questi termini perchè verranno comunque eliminati
			// M_Sys->setNullColumn(dof_IFace1[k]);  // non serve mettere a zero questi termini perchè verranno comunque eliminati
	        }

    	
}

int Problem::index_inside(size_type k, std::vector<long unsigned int>& vec) 
{

 int boolean=0;
 for (size_type i=0;i<vec.size();i++)
  {
   if (vec[i]==k)
     boolean=1;
  }

return boolean;
}








//assemblo il RHS
void Problem::assembleRHS(LinearSystem* sys)
{

	scalarVectorPtr_Type  source1;
	source1.reset(new scalarVector_Type (M_nbDOF1));
	gmm::clear(*source1);
	bulkLoad(source1, M_Bulk1, M_uFEM1, M_CoeffFEM1, M_intMethod1);
          

	scalarVectorPtr_Type  source2;
	source2.reset(new scalarVector_Type (M_nbDOF2));
	gmm::clear(*source2);
	bulkLoad(source2, M_Bulk2, M_uFEM2, M_CoeffFEM2, M_intMethod2);


// sommo i valori del rhs relativi alle funzioni definite sull'interfaccia di Omega1 a quelli relativi a Omega2
	for (size_type k=0; k<dof_IFace1.size(); k++){
		scalar_type newVal = source1->at(dof_IFace1[k]) + source2->at(dof_IFace2[k]);
		//std::cout<< "RHS riga : "<< dof_IFace2[k] + M_nbDOF1 <<" old val: "<<(M_Sys->getRHS())->at(dof_IFace2[k]+M_nbDOF1);
		source2->at(dof_IFace2[k]) = newVal;
		//M_Sys->setRHSValue(dof_IFace2[k] + M_nbDOF1, newVal);
		//std::cout<<" new val: "<<(M_Sys->getRHS())->at(dof_IFace2[k]+M_nbDOF1)<<std::endl;
		
	}
     

 //aggiungo i termini che tengono conto del salto q0 nel rhs 

scalarVectorPtr_Type q01(new scalarVector_Type(M_uFEM1.nb_dof()));
jump(q01, M_Bulk1, M_uFEM1); 
// qui suppongo di sapere q0 per vedere se il sistema funziona, poi lo prenderò dal datafile
     /* for (size_type k=0; k<M_nbDOFIFace;k++)
	{ 
            q01->at(dof_IFace1[k])=-5;  
	}*/


    for (size_type i=0; i<M_nbDOF1; i++) 
         {
		scalar_type value1=0;
	 for (size_type j=0; j< M_nbDOFIFace; j++) 
	   {
         if (index_inside(i,dof_IFace1)==0)
		{        
			value1+= (*M_Sys->getMatrix())(i,dof_IFace2[j]+M_nbDOF1)*q01->at(dof_IFace1[j]);   		   	 
		}
	   
	    }
	 source1->at(i)-= value1/2;
       } 


   // valuto q0 in tutti i nodi in modo che la numerazione per l'interfaccia rimanga uguale
        
scalarVectorPtr_Type q02(new scalarVector_Type(M_uFEM2.nb_dof()));
jump(q02, M_Bulk2, M_uFEM2); 
 /*  for (size_type k=0; k<M_nbDOFIFace;k++)
	{
          q02->at(dof_IFace2[k])=-5;	 
	} */

   for (size_type i=0; i<M_nbDOF2; i++) 
         {
         scalar_type value2=0;
	 for (size_type j=0; j< M_nbDOFIFace; j++) 
	   {
         if (index_inside(i,dof_IFace2)==0 )
		{
		  value2 += (*M_Sys->getMatrix())(i+M_nbDOF1,dof_IFace2[j]+M_nbDOF1)*q02->at(dof_IFace2[j]);
		   	 
		}
	   
	 
	   }
	 source2->at(i)+= value2/2;
       } 

	M_Sys->addSubVector(source1,0);
	M_Sys->addSubVector(source2,M_nbDOF1);



	
	/*std::cout<<" il rhs globale prima di ogni modifica è "<<std::endl;
	for (int t = 0; t<M_Sys->getRHS()->size(); t++){
		std::cout<<M_Sys->getRHS()->at(t)<<"\t";
	}*/


	
       

    
	scalarVectorPtr_Type  BCvec1;
	BCvec1.reset(new scalarVector_Type (M_nbDOF1));
	scalarVectorPtr_Type  BCvec2;
	BCvec2.reset(new scalarVector_Type (M_nbDOF2));

	stressRHS( BCvec1, M_Bulk1, &M_BC1, M_uFEM1, M_uFEM1, M_intMethod1);
	stressRHS( BCvec2, M_Bulk2, &M_BC2, M_uFEM2, M_uFEM2, M_intMethod2);


	// trovo i nodi sull'unica frontiera di Neumann che metto
	M_Sys->addSubVector(BCvec1, 0);
	M_Sys->addSubVector(BCvec2, M_nbDOF1); 


}


// risolve il sistema e estrae la soluzione
void Problem::solve()
{

       	M_Sys->solve();
      


       	M_uSol.reset(new scalarVector_Type (M_nbTotDOF));

		gmm::clear(*M_uSol);

        	M_Sys->extractSubVector(M_uSol, 0, "sol");
				//std::cout << "Size of the solution in Problem: "<< M_uSol->size() << std::endl;
				


// aggiungo di nuovo i nodi per l'interfaccia sul dominio 1 e modifico le soluzioni dove c'è la media con il valore vero anche sul dominio 2  


scalarVectorPtr_Type q01(new scalarVector_Type(M_uFEM1.nb_dof()));
jump(q01, M_Bulk1, M_uFEM1); 

  for (size_type j=0; j<M_nbDOFIFace; j++)
     {   size_type idx= dof_IFace1[j];
  	scalar_type value= M_uSol->at(dof_IFace2[j]+M_nbDOF1-M_nbDOFIFace+j)+q01->at(dof_IFace1[j])/2;
       auto it= M_uSol->insert(M_uSol->begin()+idx, value);
      // std::cout<< "la nuova dimensione di M_uSol è:" <<M_uSol->size()<<std::endl; 
    //std::cout<<"Soluzione nel dof"<<idx<<"è:"<<M_uSol->at(idx)<<std::endl;
         
	   }	




scalarVectorPtr_Type q02(new scalarVector_Type(M_uFEM2.nb_dof()));
jump(q02, M_Bulk2, M_uFEM2); 

 for (size_type j=0; j< M_nbDOFIFace; j++) 
	{
          M_uSol->at(dof_IFace2[j]+M_nbDOF1)-=q02->at(dof_IFace2[j])/2;
        }




}




// Provo a modificarla per estrarre le due soluzioni su Omega1 e Omega2 separatamente e poi anche per tenere quella globale. Definisco le due variabili M_uSol1 e M_uSol2 in questa classe Problem ed estraggo le parti da quella globale (poi eventualmente avrò bisogno di un'altro "layer" di funzioni per restituirle in output se serve)
void Problem::extractSol(scalarVectorPtr_Type destination, std::string variable)
{


	if (variable == "u1"){
			/*std::cout << "n DOFS omega1 = " << M_nbDOF1 << std::endl;
			std::cout << "sub Interval first: "<< test1.first() << "; sub Interval Last: " << test1.last() << std::endl;*/
			//destination.reset(new scalarVector_Type (M_nbDOF1));
			gmm::clear(*destination);
			gmm::copy ( gmm::sub_vector (*M_uSol, gmm::sub_interval (0,  M_nbDOF1 )), *destination); //gmm::copy(source, destination)
			std::cout << "INSIDE THE EXTRACT FUNCTION : values of M_uSol1 after "<< std::endl;
			for (size_type i=0; i < dof_IFace1.size(); i++){
				size_type idx=dof_IFace1[i];
				std::cout << M_uSol1 -> at(idx)<< "\t";
		}
	}
	else if (variable == "u2"){
		/*std::cout << "n DOFS omega2 = " << M_nbDOF2 << std::endl;
		std::cout << "n DOFS totali = " << M_nbTotDOF << std::endl;
		std::cout << "sub Interval first: "<< test2.first() << "; sub Interval Last: " << test2.last() << std::endl;*/
		//destination.reset(new scalarVector_Type (M_nbDOF2));
		gmm::clear(*destination);
		gmm::copy ( gmm::sub_vector (*M_uSol, gmm::sub_interval (M_nbDOF1,  M_nbDOF2 )), *destination); //gmm::copy(source, destination)
		std::cout << "INSIDE THE EXTRACT FUNCTION : values of M_uSol2 after "<< std::endl;
		for (size_type i=0; i < dof_IFace2.size(); i++){
				size_type idx=dof_IFace2[i];
				std::cout << M_uSol2 -> at(idx)<< "\t";
		}
	}
	else {
			destination.reset(new scalarVector_Type (M_nbTotDOF));
			gmm::clear(*destination);
			gmm::copy ( gmm::sub_vector (*M_uSol, gmm::sub_interval (0,  M_nbTotDOF )), *destination);
	}
}


// We keep the two solutions always separated since we don't have the global FEM space for the whole solution (I guess)
void Problem::exportVtk(std::string folder, std::string what)
{
	std::cout << "export "+ what <<std::endl;
	getfem::vtk_export exp(folder + "Solution_" + what + ".vtk" );
	if (what == "u1"){
		exp.exporting( *(M_uFEM1.getFEM()));
		extractSol(M_uSol1, what);
		std::cout<< "Solution extracted! "<< std::endl;
		exp.write_mesh();
		exp.write_point_data( *(M_uFEM1.getFEM()), *M_uSol1, what);
	}
	if (what == "u2"){
		exp.exporting( *(M_uFEM2.getFEM()));
		extractSol(M_uSol2, what);
		std::cout<< "Solution extracted! "<< std::endl;
		exp.write_mesh();
		exp.write_point_data( *(M_uFEM2.getFEM()), *M_uSol2, what);
	}

}



// Sets the DIRICHLET boundary conditions for the two domains in the global matrix contained in Sys. Actually it sets the rows of 0 and 1 also for the interface in the left (1) domain, since we defined the boundary condition in the Datafile as a Dirichlet one. Then we must modify these lines by adding the -1 in the columns that correspond to the basis functions on the right (2) domain.
void Problem::enforceStrongBC(size_type const domainIdx)
{
	// Set local variables according to the current domain
	//(vorrei fare delle references per non sprecare memoria, ma non posso perchè dovrei inizializzarle e non so come, quindi uso dei pointers
	FEM * M_uFEMcur;
	BC * M_BCcur;
	sizeVector_Type  M_rowsStrongBCcur;
	sizeVector_Type  M_rowsStrongBCFlagscur;

	if (domainIdx == 1){
		M_uFEMcur = &M_uFEM1;
		M_BCcur = &M_BC1;
		M_rowsStrongBCcur = M_rowsStrongBC1;
		M_rowsStrongBCFlagscur = M_rowsStrongBCFlags1;
	}
	else if (domainIdx == 2){
		M_uFEMcur = &M_uFEM2;
		M_BCcur = &M_BC2;
		M_rowsStrongBCcur = M_rowsStrongBC2;
		M_rowsStrongBCFlagscur = M_rowsStrongBCFlags2;
	}

	for ( size_type bndID = 0; bndID < M_BCcur->getDiriBD().size(); bndID++ )
	{
		dal::bit_vector quali = M_uFEMcur->getFEM()->dof_on_region( M_BCcur->getDiriBD()[bndID]);
		size_type counter=0;
		for(dal::bv_visitor i(quali); !i.finished(); ++i)
		{
			// Sono vettori di size_type cioè di indici: contengono gli indici dei dof interessati dalle BC e gli indici delle frontiere con Dirichlet, credo
			
			M_rowsStrongBCcur.push_back(i); // aggiungo gli indici dei dofs, che però non sono necessariamente in ordine!
			M_rowsStrongBCFlagscur.push_back(bndID);
			counter++;
		}
	}
                		
		scalarVectorPtr_Type q02(new scalarVector_Type(M_uFEM2.nb_dof()));
		jump(q02, M_Bulk2, M_uFEM2); 

		// Loop sui dofs
		//std::cout<<"Il valore delle BC sul dominio "<<domainIdx<<" è: "<<std::endl;
		//std::cout<< "indici dove vengo imposte bc di Dirichlet"<<std::endl;
		for (int i=0;i<M_rowsStrongBCcur.size(); ++i)
		{
			size_type ii;

			//Indice "locale" del dof
			ii = M_rowsStrongBCcur[i] ;
			bgeot::base_node where;
			where=M_uFEMcur->getFEM()->point_of_basic_dof(ii);
			//std::cout<< ii<<"\t"<<" per il nodo fisico:"<<where[0]<<","<<where[1]<<std::endl;


		// Cambia la matrice direttamente! Qui abbiamo bisogno di usare l'indice "globale" del dof : se consideriamo il secondo dominio, assumendo che la matrice in Sys sia ordinata con tutti i dof del primo dominio all'inizio e tutti gli altri alla fine, per recuperare gli indici globali delle funzioni interessate dalle BC devo aggiungere il numero di righe occupate dal primo blocco relativo al dominio1.
				 
                                M_Sys->setNullRow(ii  + (domainIdx==1 ? 0 : M_nbDOF1 ));
				M_Sys->setMatrixValue(ii  + (domainIdx==1 ? 0 : M_nbDOF1 ), ii  + (domainIdx==1 ? 0 : M_nbDOF1 ), 1);
				

				scalar_type value= M_BCcur->BCDiri(where, M_BCcur->getDiriBD()[M_rowsStrongBCFlagscur[i]]);
                                  
			   if(domainIdx==1)
				{
				M_Sys->setRHSValue(ii, value); }
                            else {
				if (index_inside(ii,dof_IFace2)==0) {
				M_Sys->setRHSValue(ii + M_nbDOF1, value); 
				} 
                                else {
				//std::cout<<"ii è:"<<ii<<std::endl;						
			        M_Sys->setRHSValue(ii  + M_nbDOF1, value+q02->at(ii)/2 ); 
				}
				
				//std::cout<<" ii è:" <<ii<<std::endl;
				//std::cout<<" i è:" <<i<<std::endl;
				//std::cout<<"punto: ("<<where[0]<<", "<<where[1]<<")\tvalore: "<<value<<std::endl; 
			}
		}


}

// To be called AFTER the enforceStrongBC --> eliminates the nodes associated to interface1 both for the matrix and the RHS
void Problem::enforceInterfaceJump(){

      M_Sys->eliminateRowsColumns(M_rowsIFace1);      // ora il sistema ha dimensione ridotta M_nbDOF1+M_nbDOF2-M_nbIFace
      M_nbTotDOF= M_nbDOF1+M_nbDOF2-M_nbDOFIFace;  

  }





// Computes the L2 and H1 errors with respect to the exact solution and saves them in the associated variables.
// NB: per prima cosa estraggo le due sottosoluzioni, anche se magari è già stato fatto dall'exportVtk.
void Problem::computeErrors(){

	extractSol(M_uSol1, "u1");
	extractSol(M_uSol2, "u2");

	scalarVectorPtr_Type DIFF1(new scalarVector_Type(M_uFEM1.nb_dof()));
	scalarVectorPtr_Type DIFF2(new scalarVector_Type(M_uFEM2.nb_dof()));

	exactSolution(DIFF1, M_Bulk1, M_uFEM1); // calcolo la soluzione esatta
	exactSolution(DIFF2, M_Bulk2, M_uFEM2);

	// Compute difference with the numerical sol
	for (size_type i=0; i<M_uFEM1.nb_dof(); i++){
		 DIFF1->at(i) -= M_uSol1->at(i);
	}

	for (size_type i=0; i<M_uFEM2.nb_dof(); i++){
		 DIFF2->at(i) -= M_uSol2->at(i);
	}

	errL2sx = getfem::asm_L2_norm(M_intMethod1, *(M_uFEM1.getFEM()), *DIFF1);
	errH1sx = getfem::asm_H1_norm(M_intMethod1, *(M_uFEM1.getFEM()), *DIFF1);
	errL2dx = getfem::asm_L2_norm(M_intMethod2, *(M_uFEM2.getFEM()), *DIFF2);
	errH1dx = getfem::asm_H1_norm(M_intMethod2, *(M_uFEM2.getFEM()), *DIFF2);

	errL2 = errL2sx + errL2dx;
	errH1 = errH1sx + errH1dx;
}



// DEBUG : print the values of the two solutions on the interface, associating the ones with the same coordinates
void Problem::printInterfaceValues(){

	std::vector<base_node> coordsSx = M_uFEM1.getDOFpoints();
	std::vector<base_node> coordsDx = M_uFEM2.getDOFpoints();
	std::cout<< " I valori della soluzione sull'interfaccia di SINISTRA e di DESTRA sono:" <<std::endl;

	for (size_type k=0; k<dof_IFace1.size(); k++){
		for (size_type h=0; h<dof_IFace2.size(); h++){

			if (gmm::vect_norm2(coordsSx[dof_IFace1[k]] - coordsDx[dof_IFace2[h]])<1.0e-9){
				std::cout<<"Coords SX: ("<<coordsSx[dof_IFace1[k]][0]<<" , "<<coordsSx[dof_IFace1[k]][1]<<")\tCoords DX: ("<<coordsDx[dof_IFace2[h]][0]<<" , "<<coordsDx[dof_IFace2[h]][1]<<")"<<std::endl;
				std::cout<<"Sx: "<< M_uSol1->at(dof_IFace1[k])<<" Dx: "<<M_uSol2->at(dof_IFace2[h])<<std::endl;
			}
		}
	}

}

// DEBUG : prints the coordinates of the nodes on the interface following the order of the vectors (dof_IFace1, dof_IFace2) that contain them, to check that they are the same (so we don't need to check every time if the coordinates are the same, but just take the element in the same position)
void Problem::printCoordinatesInterfaceNodes(){

	for (size_type k=0; k<M_nbDOFIFace; k++){
			bgeot::base_node whereSX = M_uFEM1.getFEM()->point_of_basic_dof(dof_IFace1[k]);
			bgeot::base_node whereDX = M_uFEM2.getFEM()->point_of_basic_dof(dof_IFace2[k]);

			std::cout<<" Coordinate dofs interfaccia SINISTRA: "<<whereSX[0]<<" "<<whereSX[1]<<" DESTRA: "<<whereDX[0]<< " "<<whereDX[1]<<std::endl;
	}

}
