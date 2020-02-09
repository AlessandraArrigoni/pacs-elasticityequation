#include "../include/Method2.h"

Method2::Method2(GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys):
Problem(dataFile, problem, bulk1, bulk2, dim, extSys),
M_jump1(dataFile, "bulkData/", problem, "1", "qzero"),
M_jump2(dataFile, "bulkData/", problem, "2", "qzero") {}

void Method2::assembleRHS()
{

	scalarVectorPtr_Type  source1 = std::make_shared<scalarVector_Type> (M_nbDOF1);
        bulkLoad(source1, M_uFEM1, M_uFEM1, M_source1, M_intMethod1);

	scalarVectorPtr_Type  source2 = std::make_shared<scalarVector_Type> (M_nbDOF2);
	bulkLoad(source2, M_uFEM2, M_uFEM2, M_source2, M_intMethod2);



	// Come per la matrice, "sposto" (sommandoli) anche i valori del rhs relativi alle funzioni definite sull'interfaccia di Omega1, perchè poi questi valori vengono persi quando impongo la condizione q0 con enforceInterfaceJump.
	for (size_type k=0; k<dof_IFace1.size(); k++){
		scalar_type newVal = source1->at(dof_IFace1[k]) + source2->at(dof_IFace2[k]);
		//std::cout<< "RHS riga : "<< dof_IFace2[k] + M_nbDOF1 <<" old val: "<<(M_Sys->getRHS())->at(dof_IFace2[k]+M_nbDOF1);
		source2->at(dof_IFace2[k]) = newVal;
		//std::cout<<" new val: "<<(M_Sys->getRHS())->at(dof_IFace2[k]+M_nbDOF1)<<std::endl;
	}

//aggiungo i termini che tengono conto del salto q0 nel rhs 

scalarVectorPtr_Type q01(new scalarVector_Type(M_uFEM1.nb_dof()));
jump(q01, M_jump1, M_uFEM1);      

    for (size_type i=0; i<M_nbDOF1; ++i) 
         {
	scalar_type value1=0;
	if (index_inside(i,dof_IFace1)==0){
	 for (size_type j=0; j< M_nbDOFIFace; ++j) 
	   	{       // std::cout<<" idx sull'interfaccia è:"<<dof_IFace1[j]<<std::endl;
			value1+= (*(M_Sys.getMatrix()))(i,dof_IFace2[j]+M_nbDOF1)*q01->at(dof_IFace1[j]);   	 
		}
	   
	  source1->at(i)-= value1/2;
       } 
}

 // valuto q0 in tutti i nodi in modo che la numerazione per l'interfaccia rimanga uguale
        
scalarVectorPtr_Type q02(new scalarVector_Type(M_uFEM2.nb_dof()));
jump(q02, M_jump2, M_uFEM2); 


   for (size_type i=0; i<M_nbDOF2; ++i) 
         {
         scalar_type value2=0;

      if (index_inside(i,dof_IFace2)==0 )
		{
	 for (size_type j=0; j< M_nbDOFIFace; ++j) 
	   { 	
		  value2 += (*(M_Sys.getMatrix()))(i+M_nbDOF1,dof_IFace2[j]+M_nbDOF1)*q02->at(dof_IFace2[j]);
		  
	   }
	 
	 source2->at(i)+= value2/2; 
	  
       } 
}

	M_Sys.addSubVector(source1,0);
	M_Sys.addSubVector(source2,M_nbDOF1);




  // To include the Neumann boundary conditions (maybe then we will need to change this part adding the jump too - or defining it in the dataFile in a proper way)

       	scalarVectorPtr_Type BCvec1 = std::make_shared<scalarVector_Type> (M_nbDOF1);
	stressRHS( BCvec1, M_Bulk1, M_BC1, M_uFEM1, M_uFEM1, M_intMethod1);

	scalarVectorPtr_Type BCvec2 = std::make_shared<scalarVector_Type> (M_nbDOF2);
	stressRHS( BCvec2, M_Bulk2, M_BC2, M_uFEM2, M_uFEM2, M_intMethod2);

	M_Sys.addSubVector(BCvec1, 0);
	M_Sys.addSubVector(BCvec2, M_nbDOF1);

}




void Method2::enforceStrongBC(size_type const domainIdx)
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

size_type Qdim = M_uFEMcur->getFEM().get_qdim();

		for ( size_type bndID = 0; bndID < M_BCcur->getDiriBD().size(); bndID++ )
		{
			dal::bit_vector quali = M_uFEMcur->getFEM().dof_on_region( M_BCcur->getDiriBD()[bndID]);
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
jump(q02, M_jump2, M_uFEM2);
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
				where = M_uFEMcur->getFEM().point_of_basic_dof(ii + j);	
                                M_Sys.setNullRow(ii + j  + (domainIdx==1 ? 0 : M_nbDOF1 ));
				M_Sys.setMatrixValue(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), ii  + j + (domainIdx==1 ? 0 : M_nbDOF1 ), 1);

				scalar_type value= M_BCcur->BCDiri(where,j,M_BCcur->getDiriBD()[M_rowsStrongBCFlagscur[i]]);
                                  
				
				 if(domainIdx==1)
				{
				M_Sys.setRHSValue(ii+j, value); }
                            else {
				if (index_inside(ii+j,dof_IFace2)==0) {
				M_Sys.setRHSValue(ii + j +M_nbDOF1, value); 
				} 
                                else {
				//std::cout<<"ii è:"<<ii<<std::endl;						
			        M_Sys.setRHSValue(ii  + j +M_nbDOF1, value+q02->at(ii+j)/2 ); 
				}
			       }
			 } // fine loop su j (2 dofs associati a un punto fisico)
										
		}								
										 
}


void Method2::treatIFaceDofs(){

      M_Sys.eliminateRowsColumns(dof_IFace1);      // ora il sistema ha dimensione ridotta M_nbDOF1+M_nbDOF2-M_nbIFace
      M_nbTotDOF= M_nbDOF1+M_nbDOF2-M_nbDOFIFace;  

  }


size_type Method2::index_inside(size_type k, std::vector<long unsigned int>& vec) 
{

 size_type boolean = 0;
 for (size_type i=0;i<vec.size();i++)
  {
   if (vec[i]==k)	

     boolean=1;
  }

return boolean;
}



void Method2::solve()
{

  M_uSol.resize(M_nbTotDOF, 0.0);
  M_uSol1.resize(M_nbDOF1, 0.0);
  M_uSol2.resize(M_nbDOF2, 0.0);

       	M_Sys.solve();

       	M_Sys.extractSubVector(M_uSol, 0, "sol");
				std::cout << "Size of the solution in Problem: "<< M_uSol.size() << std::endl;
				// The content of the variable M_sol contained in the LinearSystem M_Sys (actually it is all a matter of pointers...) is copied into the variable M_uSol of our class Problem.

// aggiungo di nuovo i nodi per l'interfaccia sul dominio 1 e modifico le soluzioni dove c'è la media con il valore vero anche sul dominio 2  


scalarVectorPtr_Type q01(new scalarVector_Type(M_uFEM1.nb_dof()));
jump(q01, M_jump1, M_uFEM1); 


  for (size_type j=0; j<M_nbDOFIFace; j++)
	
     {   size_type idx= dof_IFace1[j];
  	scalar_type value= M_uSol.at(dof_IFace2[j]+M_nbDOF1-M_nbDOFIFace+j)+q01->at(idx)/2; // valore della media trovato + il salto/2
       auto it= M_uSol.insert(M_uSol.begin()+idx, value);
      // std::cout<< "la nuova dimensione di M_uSol è:" <<M_uSol->size()<<std::endl; 
    //std::cout<<"Soluzione nel dof"<<idx<<"è:"<<M_uSol->at(idx)<<std::endl;
             
	   }	

scalarVectorPtr_Type q02(new scalarVector_Type(M_uFEM2.nb_dof()));
jump(q02, M_jump2, M_uFEM2); 

 for (size_type j=0; j< M_nbDOFIFace; j++) 
	
	{
          M_uSol.at(dof_IFace2[j]+M_nbDOF1)-=q02->at(dof_IFace2[j])/2; // il valore della media - il salto/2
        }


  
}
