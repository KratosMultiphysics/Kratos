#if !defined(KRATOS_SPLIT_TRIANGLE)
#define  KRATOS_SPLIT_TRIANGLE

//input parametrs:
//edges --> int c array of size 
//t     --> int c array of size 12 (3*4)
//nel = number of elements in the subdivision
//output parameters: true->splitting needed    false-->no splitting needed
//nint Nodo interno dentro del triangulo
         
bool Split_Triangle(const int* edges, int* t, int* nel, int* splitted_edges, int* nint) 
          {     
            *splitted_edges = 0;     
	     bool topology[3];
	     topology[0] = false; 
	     topology[1] = false;
	     topology[2] = false;
	    for (unsigned int i = 0; i < 3; i++)
	        {
	           if (edges[i] > 2)
		   {  topology[i] = true;
	              *splitted_edges = *splitted_edges + 1;
		   }
		 }
	    	  
	    if (*splitted_edges == 0 && *nint == 0) 
	    {
	    //no splitting needed
	    *nel = 1;
	    t[0] = 0;
	    t[1] = 1;
	    t[2] = 2;
	    return false;
	    }
	    
	   ///WARNING = Caso de un nodo central dentro del elemento 
	   else if (*splitted_edges == 0 && *nint == 1) 
	    {
	    *nel = 3;
	   
	    t[0] = 3;
	    t[1] = 0;
	    t[2] = 1;
	   
	    t[3] = 3;
	    t[4] = 1;
	    t[5] = 2;
	   
	    t[5] = 3;
	    t[5] = 2;
	    t[5] = 0;
	    return true;
	    }
                
                     
            else if(*splitted_edges==1)
               {      
		 *nel = 2;
                 /// caso 1
                 if(topology[0]== true)
                  {
                   t[0] = 3;
		   t[1] = 2;
		   t[2] = 0; 
		   
                   t[3] = 3; 
		   t[4] = 1;
		   t[5] = 2;
                  } 

                /// caso 2
              else if( topology[1]== true)
                  {
                   t[0] = 4; 
		   t[1] = 0; 
		   t[2] = 1; 
		   
                   t[3] = 4;
		   t[4] = 2; 
		   t[5] = 0;
                  }

                /// caso 3 
              else if( topology[2]== true)
                    {
                      t[0] = 5;  
		      t[1] = 1;
		      t[2] = 2; 
		      
                      t[3] = 5;  
		      t[4] = 0; 
		      t[5] = 1;
                    }

                return true;                    

               }
             
             
             else if(*splitted_edges==2)
               {
		     *nel = 3;
                     /// caso 4  
                     if( topology[0]== true &&  topology[1]== true )
                       {
			   if(edges[2]==0) // si colapso al nodo 0 local
			   {
                              t[0] =  4;      
			      t[1] =  3;  
			      t[2] =  1; 
			      
                              t[3] =  4;      
			      t[4] =  0; 
			      t[5] =  3;
			      
                              t[6] =  4; 
			      t[7] =  2;    
			      t[8] =  0;                                                                  
                           }
                           
                           else if(edges[2]==2) // si colapso al nodo 2 local
			   {
			      t[0] =  4;      
			      t[1] =  3;  
			      t[2] =  1; 
			      
                              t[3] =  4;      
			      t[4] =  2; 
			      t[5] =  3;
			      
                              t[6] =  3; 
			      t[7] =  2;    
			      t[8] =  0;  
			   }
			     
		      }
	              
                     /// caso 5 
                     else if(topology[1]== true &&  topology[2]== true)
		     {
			   if(edges[0]==0) // si colapso al nodo 0 local
			   {
                              t[0] =  5;      
			      t[1] =  4;  
			      t[2] =  2; 
			      
                              t[3] =  5;      
			      t[4] =  0; 
			      t[5] =  4;
			      
                              t[6] =  4; 
			      t[7] =  0;    
			      t[8] =  1;                                                                  
                           }
                           
                           else if(edges[0]==1) // si colapso al nodo 2 local
			   {
			      t[0] =  5;      
			      t[1] =  4;  
			      t[2] =  2; 
			      
                              t[3] =  5;      
			      t[4] =  1; 
			      t[5] =  4;
			      
                              t[6] =  5; 
			      t[7] =  0;    
			      t[8] =  1;  
			   }                         
                       } 
                       
                     
                     /// caso 3 
                     else if(topology[0]== true &&  topology[2]== true)
                       {          
                          if(edges[1]==1) // si colapso al nodo 0 local
			   {
                              t[0] =  5;      
			      t[1] =  0;  
			      t[2] =  3; 
			      
                              t[3] =  5;      
			      t[4] =  3; 
			      t[5] =  1;
			      
                              t[6] =  5; 
			      t[7] =  1;    
			      t[8] =  2;                                                                  
                           }
                           
                           else if(edges[1]==2) // si colapso al nodo 2 local
			   {
			      t[0] =  5;      
			      t[1] =  0;  
			      t[2] =  3; 
			      
                              t[3] =  5;      
			      t[4] =  3; 
			      t[5] =  2;
			      
                              t[6] =  3; 
			      t[7] =  1;    
			      t[8] =  2;  
			   }  	 
			 
                       }
                       
                   return true;  
               }
 
                     else if(*splitted_edges==3)
                        {
                              *nel =  4;
			      t[0] =  5;      
			      t[1] =  0;  
			      t[2] =  3; 

			      t[3] =  5;      
			      t[4] =  3; 
			      t[5] =  4;

			      t[6] =  4; 
			      t[7] =  3;    
			      t[8] =  1;
			      
			      t[9]  =  5; 
			      t[10] =  4;    
			      t[11] =  2;
			      
			      return true; 
		      }
		     else 
                      {      
                              return false;  
                      }
                
	}

#endif // KRATOS_SPLIT_TRIANGLE  defined 

