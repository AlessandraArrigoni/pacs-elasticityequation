#include "../include/Problem.h"

Problem::Problem (const GetPot& dataFile1,const GetPot& dataFile2, Bulk* bulk1, Bulk* bulk2):
			M_Bulk1(bulk1), M_Bulk2(bulk2),
			M_BC1(dataFile1, "elasticity/", "bulkData/"), // change the dataFile accordingly and set it
			M_BC2(dataFile2, "elasticity/", "bulkData/"), // change the dataFile accordingly and set it
			interfaceIdx1(1),// change the dataFile accordingly and set it: devo capire come leggere il valore direttamente dal file!
			interfaceIdx2(3), // adesso le ho messe brutali così, poi andranno lette da file o altro.
			M_uFEM1( bulk1->getMesh(), dataFile1, "elasticity/", "Sol", "bulkData/", 2),// the last parameter is the dimesion of the space (2 = vector)
			M_uFEM2( bulk2->getMesh(), dataFile2, "elasticity/", "Sol", "bulkData/", 2),// the last parameter is the dimesion of the space (2 = vector)
			M_CoeffFEM1( bulk1->getMesh(), dataFile1, "elasticity/", "Coeff", "bulkData/"), // change the dataFile accordingly and set it
			M_CoeffFEM2( bulk2->getMesh(), dataFile2, "elasticity/", "Coeff", "bulkData/"), // change the dataFile accordingly and set it
			M_intMethod1(*(bulk1->getMesh()) ),
			M_intMethod2(*(bulk2->getMesh()) ),
		M_Sys()
{

   M_nbDOF1 = M_uFEM1.nb_dof();
	 M_nbDOF2 = M_uFEM2.nb_dof();
	 M_nbTotDOF = M_nbDOF1 + M_nbDOF2;
	 std::cout<<"Il numero totale di dofs è "<<M_nbTotDOF<<std::endl;

	 Qdim = M_uFEM1.getFEM()->get_qdim(); // Dimensione dello spazio FEM (1 = scalare, 2 = vettoriale)
	 std::cout<<"Le dimensioni degli spazi FEM sono: "<<Qdim<<" per Omega1 e "<< M_uFEM2.getFEM()->get_qdim()<<" per Omega2 e "<<M_CoeffFEM1.getFEM()->get_qdim()<<" per i coefficienti"<<std::endl;

	 // Per ora faccio tutto sdoppiato, anche se poi probabilmente il metodo di integrazione sarà lo stesso sui 2 domini. Per ora ho solo copiato e incollato.
   std::string intMethod1(dataFile1 ( std::string("bulkData/elasticity/integrationMethod" ).data (), "IM_TRIANGLE(6)" ) );
	 std::string intMethod2(dataFile2 ( std::string("bulkData/elasticity/integrationMethod" ).data (), "IM_TRIANGLE(6)" ) );

   M_intMethod1.set_integration_method(bulk1->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod1) );
	 M_intMethod2.set_integration_method(bulk2->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod2) );

   M_Bulk1->getData()->setLambda(M_CoeffFEM1.getDOFpoints());  //valuta i coefficienti di Lamè negli elementi della mesh
	 M_Bulk1->getData()->setMu(M_CoeffFEM1.getDOFpoints());
	 M_Bulk2->getData()->setLambda(M_CoeffFEM2.getDOFpoints());  //valuta i coefficienti di Lamè negli elementi della mesh
	 M_Bulk2->getData()->setMu(M_CoeffFEM2.getDOFpoints());

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
 	{		 M_rowsIFace1.push_back(i);  }

	fromBitVectorToStdVector ( dal_dof_IFace1, dof_IFace1);
	/*std::cout<<"Gli indici sull'interfaccia di SINISTRA sono: "<<std::endl;
	for (size_type j=0; j<dof_IFace1.size(); j++){
	 std::cout<< dof_IFace1[j]<<"  ";
 }*/
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
 	{		M_rowsIFace2.push_back(i);	 }

	fromBitVectorToStdVector ( dal_dof_IFace2, dof_IFace2 );
	/*std::cout<<"Gli indici sull'interfaccia di DESTRA sono: "<<std::endl;
	for (size_type j=0; j<dof_IFace2.size(); j++){
	 std::cout<< dof_IFace2[j]<<"  ";
 }*/

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
	linearElasticity(A1, M_Bulk1, M_uFEM1, M_CoeffFEM1, M_intMethod1);

	A2.reset(new sparseMatrix_Type (M_nbDOF2, M_nbDOF2));
	gmm::clear(*A2);
	linearElasticity(A2, M_Bulk2, M_uFEM2, M_CoeffFEM2, M_intMethod2);

	M_Sys->addSubMatrix(A1, 0, 0);
	M_Sys->addSubMatrix(A2, M_nbDOF1, M_nbDOF1);

	// Get all the rows associated to the dofs on the interface;
	for (size_type k=0; k<dof_IFace1.size(); k++){
		sparseMatrixPtr_Type curRow(new sparseMatrix_Type(1,M_nbTotDOF));

		M_Sys->extractSubMatrix( dof_IFace1[k], 1, 0, M_nbTotDOF, curRow );
		//std::cout<<"ho estratto la riga numero "<<dof_IFace1[k] <<std::endl;
		M_Sys->addSubMatrix(curRow, dof_IFace2[k] + M_nbDOF1, 0);
		//std::cout<<"ho aggiunto la riga numero "<<dof_IFace2[k]  + M_nbDOF1<<std::endl;

		/*
		if(k==3){
			// DEBUG : stampa la 4 (a caso) riga estratta e controllo che poi sia sommata a quella corrispondente nel sistema globale
			std::cout<<"La quarta riga estratta dalla matrice di sx è "<<std::endl;
			for(size_type j=0; j<M_nbTotDOF; j++){
				std::cout<< (*curRow)(0,j)<<"  ";
			}
			std::cout<< "\nLa quarta riga cambiata della matrice globale è "<<std::endl;
			for(size_type j=0; j<M_nbTotDOF; j++){
				std::cout<<(*(M_Sys->getMatrix()))(dof_IFace2[3] + M_nbDOF1,j)<<"  ";
			}
		}*/


	}


}

//assemblo il RHS
void Problem::assembleRHS(LinearSystem* sys)
{

	scalarVectorPtr_Type  source1;
	source1.reset(new scalarVector_Type (M_nbDOF1));
	gmm::clear(*source1);
	bulkLoad(source1, M_Bulk1, M_uFEM1, M_uFEM1, M_intMethod1);

	scalarVectorPtr_Type  source2;
	source2.reset(new scalarVector_Type (M_nbDOF2));
	gmm::clear(*source2);
	bulkLoad(source2, M_Bulk2, M_uFEM2, M_uFEM2, M_intMethod2);

	M_Sys->addSubVector(source1,0);
	M_Sys->addSubVector(source2,M_nbDOF1);


	// Come per la matrice, "sposto" (sommandoli) anche i valori del rhs relativi alle funzioni definite sull'interfaccia di Omega1, perchè poi questi valori vengono persi quando impongo la condizione q0 con enforceInterfaceJump.
	for (size_type k=0; k<dof_IFace1.size(); k++){
		scalar_type newVal = (M_Sys->getRHS())->at(dof_IFace1[k]) + (M_Sys->getRHS())->at(dof_IFace2[k] + M_nbDOF1);
		//std::cout<< "RHS riga : "<< dof_IFace2[k] + M_nbDOF1 <<" old val: "<<(M_Sys->getRHS())->at(dof_IFace2[k]+M_nbDOF1);
		M_Sys->setRHSValue(dof_IFace2[k] + M_nbDOF1, newVal);
		//std::cout<<" new val: "<<(M_Sys->getRHS())->at(dof_IFace2[k]+M_nbDOF1)<<std::endl;
	}


  // To include the Neumann boundary conditions (maybe then we will need to change this part adding the jump too - or defining it in the dataFile in a proper way)
	scalarVectorPtr_Type  BCvec1;
	BCvec1.reset(new scalarVector_Type (M_nbDOF1));
	scalarVectorPtr_Type  BCvec2;
	BCvec2.reset(new scalarVector_Type (M_nbDOF2));

	stressRHS( BCvec1, M_Bulk1, &M_BC1, M_uFEM1, M_uFEM1, M_intMethod1);
	stressRHS( BCvec2, M_Bulk2, &M_BC2, M_uFEM2, M_uFEM2, M_intMethod2);

	// trovo i nodi sull'unica frontiera di Neumann che metto
	M_Sys->addSubVector(BCvec1, 0);
	M_Sys->addSubVector(BCvec2, M_nbDOF1);

	/*
	std::cout<<"\nIl rhs globale dopo l'imposizione della BC di Neumann è "<<std::endl;
	for (int t = 0; t<M_Sys->getRHS()->size(); t++){
		std::cout<<M_Sys->getRHS()->at(t)<<"\t";
	}*/
}


// risolve il sistema e estrae la soluzione
void Problem::solve()
{

       	M_Sys->solve();

       	M_uSol.reset(new scalarVector_Type (M_nbTotDOF));

				gmm::clear(*M_uSol);

       	M_Sys->extractSubVector(M_uSol, 0, "sol");
				std::cout << "Size of the solution in Problem: "<< M_uSol->size() << std::endl;
				// The content of the variable M_sol contained in the LinearSystem M_Sys (actually it is all a matter of pointers...) is copied into the variable M_uSol of our class Problem.
}

// Provo a modificarla per estrarre le due soluzioni su Omega1 e Omega2 separatamente e poi anche per tenere quella globale. Definisco le due variabili M_uSol1 e M_uSol2 in questa classe Problem ed estraggo le parti da quella globale (poi eventualmente avrò bisogno di un'altro "layer" di funzioni per restituirle in output se serve)
void Problem::extractSol(scalarVectorPtr_Type destination, std::string variable)
{

/*
	gmm::sub_interval test2(M_nbDOF1, M_nbDOF2 );
	gmm::sub_interval test1(0,  M_nbDOF1);
*/
	if (variable == "u1"){
			/*std::cout << "n DOFS omega1 = " << M_nbDOF1 << std::endl;
			std::cout << "sub Interval first: "<< test1.first() << "; sub Interval Last: " << test1.last() << std::endl;*/
			//destination.reset(new scalarVector_Type (M_nbDOF1));
			gmm::clear(*destination);
			gmm::copy ( gmm::sub_vector (*M_uSol, gmm::sub_interval (0,  M_nbDOF1 )), *destination); //gmm::copy(source, destination)
			//std::cout << "INSIDE THE EXTRACT FUNCTION : values of M_uSol1 after "<< std::endl;
			//for (size_type i=0; i < 5; i++){
			//	std::cout << M_uSol1 -> at(i)<< "\t";
		// }
	}
	else if (variable == "u2"){
		/*std::cout << "n DOFS omega2 = " << M_nbDOF2 << std::endl;
		std::cout << "n DOFS totali = " << M_nbTotDOF << std::endl;
		std::cout << "sub Interval first: "<< test2.first() << "; sub Interval Last: " << test2.last() << std::endl;*/
		//destination.reset(new scalarVector_Type (M_nbDOF2));
		gmm::clear(*destination);
		gmm::copy ( gmm::sub_vector (*M_uSol, gmm::sub_interval (M_nbDOF1,  M_nbDOF2 )), *destination); //gmm::copy(source, destination)
		//std::cout << "INSIDE THE EXTRACT FUNCTION : values of M_uSol2 after "<< std::endl;
		//for (size_type i=0; i < 5; i++){
		//	std::cout << M_uSol2-> at(i)<< "\t";
		//}
	}
	else {
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
			dal::bit_vector quali = M_uFEMcur->getFEM()->dof_on_region( M_BCcur->getDiriBD()[bndID]); // QUESTO PUÒ ESSERE UTILE PER OTTENERE I DOFS SULL'INTERFACCIA SE RIESCO A CAPIRE COME! getDiriBD restituisce un vettore che contiene gli indici (1,2,3,4) dei lati interessati da una BC di Dirichlet.
			size_type counter=0;
			for(dal::bv_visitor i(quali); !i.finished(); ++i)
			{
				// Sono vettori di size_type cioè di indici: contengono gli indici dei dof interessati dalle BC e gli indici delle frontiere con Dirichlet, credo
				M_rowsStrongBCcur.push_back(i); // aggiungo gli indici dei dofs, che però non sono necessariamente in ordine!
				M_rowsStrongBCFlagscur.push_back(bndID);
				counter++;
			}
		}

		// Loop sui dofs: ora questo mi serve solo per creare il vettore con i dati da passare alla funzione assembling_Dirichlet_condition: l'idea è che il vettore sia nullo ovunque, tranne nei dof sulle frontiere. NON SONO SICURA! FORSE DEVO GIÀ PASSARE QUESTO SOMMATO AL RHS CHE GIÀ HO. In realtà penso di no perchè, guardando il codice, nei dofs di rhs che non sono sulle frontiere toglie A*ug, ma se definisco ug come nulla su quei dof, dovrei essere a posto e non avere nessuna modifica al rhs.

		//IPOTESI MOLTO FORTE DA CONTROLLARE, altrimenti devo trovare il modo di fare il loop sui nodi fisici che ci sono sulle frontiere e recuperare poi i dofs (sperando che il primo sia relativo alla componente x e il secondo alla y): quando recupero i dof_on_region, mi restituisce accoppiati i due relativi a ogni nodo fisico, quindi poi anche qui quando faccio il loop su quelli sono sicura che avrò le stesse coordinate del where in due step consecutivi.

		//std::cout<<"Il valore delle BC sul dominio "<<domainIdx<<" è: "<<std::endl;
		for (size_type i = 0; i < M_rowsStrongBCcur.size(); i += Qdim)
		{
			size_type ii; // é l'indice "locale" del primo dof associato a un punto fisico (spero!)
			ii = M_rowsStrongBCcur[i] ;

										// DEBUG
										//bgeot::base_node whereX, whereY;
			for(size_type j = 0; j < Qdim; j++)
			{
				bgeot::base_node where;
				where = M_uFEMcur->getFEM()->point_of_basic_dof(ii + j);
				scalar_type value = M_BCcur->BCDiri(where, j, M_BCcur->getDiriBD()[M_rowsStrongBCFlagscur[i]]);

										// DEBUG
										//if (j==0){whereX = where;}
										//if (j==1){whereY = where;}

				M_Sys->setNullRow(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ));
				M_Sys->setMatrixValue(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), 1);

				M_Sys->setRHSValue(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), value);

			} // fine loop su j (2 dofs associati a un punto fisico)

										//DEBUG
										//std::cout<<"Nuova coppia di dofs"<<std::endl;
										//std::cout<<"punto dal dof X : ("<<whereX[0]<<", "<<whereX[1]<<")"<<std::endl;
										//std::cout<<"punto dal dof Y : ("<<whereY[0]<<", "<<whereY[1]<<")"<<std::endl;
		}

		/*
		std::cout<<"\nIl rhs globale dopo l'imposizione della BC di Dirichlet sul dominio "<<domainIdx<<" è "<<std::endl;
		for (int t = 0; t<M_Sys->getRHS()->size(); t++){
			std::cout<<M_Sys->getRHS()->at(t)<<"\t";
		}
		std::cout<<std::endl;
		*/
}

// To be called AFTER the enforceStrongBC since that function sets the whole row to 0 (even the one associated to the interface!)
// devo recuperare gli indici relativi all'interfaccia di entrambi i domini e poi modificare le righe della matrice relative a omega1: in realtà potrei solo aggiungere i -1 in corrispondenza delle colonne di omega2, visto che enforceStrongBC modifica già le righe corrispondenti, ma meglio rifare tutto.
void Problem::enforceInterfaceJump(){

	size_type idxIFaceInDiriBD;
 	for (size_type j=0; j<M_BC1.getDiriBD().size(); j++){
 		if (M_BC1.getDiriBD()[j] == interfaceIdx1) {idxIFaceInDiriBD = j;}
 	}

	// Loop on the interface dofs and change the matrix
	//std::cout<<"I valori sull'INTERFACCIA sono: "<<std::endl;
	for (size_type k = 0 ; k < M_nbDOFIFace ; k += Qdim){

		// Global indices in the two domains: vectors dof_IFace1 and dof_IFace2 are initialized by the constructor. I hope, once again, that the 2 dofs associated to the same physical node occupy two consecutive places in the vector!
		size_type idx1 = dof_IFace1[k];
		size_type idx2 = dof_IFace2[k] + M_nbDOF1;

		for (size_type j = 0; j < Qdim; j++)
		{
			// Get the coordinates of the points
			bgeot::base_node where;
			where = M_uFEM1.getFEM()->point_of_basic_dof(idx1 + j);

			// Change the matrix in M_Sys
			M_Sys->setNullRow(idx1 + j);
			M_Sys->setMatrixValue(idx1 + j, idx1 + j, 1);
			M_Sys->setMatrixValue(idx1 + j, idx2 + j, -1);

			// Change the RHS
			scalar_type value= M_BC1.BCDiri(where, j, M_BC1.getDiriBD()[idxIFaceInDiriBD]);
			M_Sys->setRHSValue(idx1 + j, value);

		} // fine loop su j (2 dofs associati a un punto fisico)

				//std::cout<<"punto: ("<<where[0]<<", "<<where[1]<<")\tvalore: "<<value<<std::endl;
	}
	/*
	std::cout<<"Le righe che ho cambiato con la condizione di interfaccia sono : "<<std::endl;
	for (size_type k = 0 ; k<M_nbDOFIFace ; k++){
		std::cout<<dof_IFace1[k]<<" ";
	}
	std::cout<<"\n\nI valori del rhs nelle righe associate all'interfaccia di omega2 sono : "<<std::endl;
	for (size_type k = 0 ; k<M_nbDOFIFace ; k++){
		std::cout<<(M_Sys->getRHS())->at(dof_IFace2[k]+M_nbDOF1)<<" ";
	}*/
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

	// Compute difference with the numerical solution
	for (size_type i=0; i<M_uFEM1.nb_dof(); i++){
		 DIFF1->at(i) -= M_uSol1->at(i);
	}

		// DEBUG
		/*std::cout<<"\nla soluzione NUMERICA su omega 1 vale : "<<std::endl;
		for (size_type i=0; i<M_uFEM1.nb_dof(); i++){
			std::cout<<M_uSol1->at(i)<<"\t";
		}*/


	for (size_type i=0; i<M_uFEM2.nb_dof(); i++){
		 DIFF2->at(i) -= M_uSol2->at(i);
	}

	//DEBUG
	std::cout<< "\nvalore variabili errore prima di calcolarlo: L2 "<<errL2<<" e H1 "<<errH1<<std::endl;

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
