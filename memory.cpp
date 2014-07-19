#include "common.h"
#include "sim.h"
#include "trace.h" 
#include "memory.h" 
#include "knob.h"
#include "all_knobs.h"

#include <stdlib.h>

extern UINT64 dram_row_buffer_hit_count;
extern UINT64 dram_row_buffer_miss_count;
/*******************************************************************/
/* Memory Related frame */ 
/*******************************************************************/
/* You do not need to change the code in this part                 */ 


void memory_c::dprint_queues(void) 
{
  // traverse the queue and print out all the values inside the queue 

  /* print out mshr entries */ 
  
  list<m_mshr_entry_s *>::const_iterator cii; 
  
  /* do not modify print MSHR print statement here. These will be used for grading */ 

  cout <<"** cycle_count: " << cycle_count << endl; 
  cout <<"***MSHR entry print outs: mshr_size: " <<  m_mshr.size() << endl; 


  for (cii= m_mshr.begin() ; cii != m_mshr.end(); cii++) {
    
    int jj = 0; 
    cout << "mshr_entry_num: " << jj;
    jj++; 

    m_mshr_entry_s* entry = (*cii);
    
    if (entry->valid) {
      mem_req_s *req = entry->m_mem_req; 
      
      cout <<" mem_req_id: " << req->m_id; 
      cout <<" mem_addr: " << req->m_addr; 
      cout <<" size: " << (int)req->m_size; 
      cout <<" insert_time: " << entry->insert_time;
      cout <<" rdy_cycle: " << req->m_rdy_cycle; 
      cout <<" dirty: " << req->m_dirty; 
      cout <<" done: " << req->m_done; 
      cout <<" state: " << req->m_state; 
      cout <<" type:  " << req->m_type; 
      cout <<" ops_size: " << entry->req_ops.size(); 
      
      if (entry->req_ops.size()) {
	list<Op *>::const_iterator cii; 
	int kk = 0; 
	for (cii = entry->req_ops.begin() ; cii != entry->req_ops.end(); cii++) {
	  
	  Op *m_op = (*cii);
	  
	  if(m_op->mem_type == MEM_LD) {
	    printf("op[%d]:LD id: %lu ", kk, (uint64_t)m_op->inst_id);
	  }
	  else 
	    printf("op[%d]:ST id: %lu", kk, (uint64_t)m_op->inst_id);
	  	  
	}
	kk++;

      }
     
      cout << endl; 
    }
  }
  
  // print queues 

  
  cout <<"***DRAM_IN_QUEUE entry print outs ****" << endl; 
  
  list<mem_req_s *>::const_iterator cqii; 
  for (cqii = dram_in_queue.begin() ; cqii != dram_in_queue.end(); cqii++) {
    mem_req_s *req = *cqii; 
    
    cout << " req_id: " << req->m_id; 
    cout << " addr: " << req->m_addr; 
    cout << " rdy_cycle: " << req->m_rdy_cycle; 
    cout <<" state: " << req->m_state; 
    cout <<" type:  " << req->m_type; 
    cout << " ||| " ; 
  }
  
  cout << endl; 
  // end printing dram in queue 
  

  // print dram scheduler queues 
  
  list<mem_req_s *>::const_iterator clii; 
  
  cout <<"***DRAM_SCH_QUEUE entry print outs ****" << endl; 

  for (int bb = 0; bb < m_dram_bank_num; bb++) { 
    cout <<"***DRAM_SCH_QUEUE BANK[" << bb << "]" << endl; 
    
    for (clii = dram_bank_sch_queue[bb].begin(); 
	 clii != dram_bank_sch_queue[bb].end(); 
	 clii++) { 
      
      mem_req_s *req = *clii; 
      
      cout << " req_id: " << req->m_id; 
      cout << " addr: " << req->m_addr; 
      cout << " rdy_cycle: " << req->m_rdy_cycle; 
      cout <<" state: " << req->m_state; 
      cout <<" type:  " << req->m_type; 
      cout << " ||| " ; 
    
    }
    cout << endl; 
  }
  cout << endl; 
  // ending print dram scheduler 

  
  // print dram out queue 
  
  cout <<"***DRAM_OUT_QUEUE entry print outs ****" << endl; 
  for (cqii= dram_out_queue.begin() ; cqii != dram_out_queue.end(); cqii++) {
    mem_req_s *req = *cqii; 
    
    cout << " req_id: " << req->m_id; 
    cout << " addr: " << req->m_addr; 
    cout << " rdy_cycle: " << req->m_rdy_cycle; 
    cout <<" state: " << req->m_state; 
    cout <<" type:  " << req->m_type; 
    cout << " ||| " ; 
      
  }
  
  cout << endl; 
  
  
}

/*******************************************************************/
/* print dram status. debugging help function 
/*******************************************************************/

void memory_c::dprint_dram_banks()
{
  
  cout << " DRAM_BANK_STATUS cycle_count "<< cycle_count << endl; 

  for (int ii = 0; ii <m_dram_bank_num; ii++) {
    printf("bank_num[%d]: row_id:%ld rdy_cycle:%lu \n",
	   ii, dram_bank_open_row[ii], (uint64_t)dram_bank_rdy_cycle[ii]); 
  }
  
}


/*******************************************************************/
/* Initialize the memory related data structures                   */
/*******************************************************************/

void memory_c::init_mem()
{
  /* For Lab #2, you do not need to modify this code */

  /* you can add other code here if you want */ 
  

  /* init mshr */ 
  m_mshr_size = KNOB(KNOB_MSHR_SIZE)->getValue(); 
  m_dram_bank_num = KNOB(KNOB_DRAM_BANK_NUM)->getValue(); 
  m_block_size = KNOB(KNOB_BLOCK_SIZE)->getValue();

  
  for (int ii = 0 ; ii < m_mshr_size; ii++) 
  {

    m_mshr_entry_s* entry = new m_mshr_entry_s; 
    entry->m_mem_req = new mem_req_s;  // create a memory rquest data structure here 
    
    entry->valid = false; 
    entry->insert_time = 0; 
    m_mshr_free_list.push_back(entry); 
  }
  
  /* init DRAM scheduler queues */ 

  dram_in_queue.clear();
  dram_out_queue.clear();
  
  dram_bank_sch_queue = new list<mem_req_s*>[m_dram_bank_num]; 
  dram_bank_open_row = new int64_t[m_dram_bank_num]; 
  for(int ii=0;ii<m_dram_bank_num;ii++)
  {
    dram_bank_open_row[ii]=-1;
  }
  dram_bank_rdy_cycle =  new uint64_t[m_dram_bank_num]; 
  for(int ii=0;ii<m_dram_bank_num;ii++)
    {
	  dram_bank_rdy_cycle[ii]=0;
    }
  
  
}

/*******************************************************************/
/* free MSHR entry (init the corresponding mshr entry  */ 
/*******************************************************************/

void memory_c::free_mshr_entry(m_mshr_entry_s *entry)
{
  entry->valid = false; 
  entry->insert_time = 0; 
  m_mshr_free_list.push_back(entry); 

}

/*******************************************************************/
/* This function is called every cycle                             */ 
/*******************************************************************/

void memory_c::run_a_cycle() 
{
	
  
  if (KNOB(KNOB_PRINT_MEM_DEBUG)->getValue()) {
    dprint_queues();
    dprint_dram_banks();
  }
  
  /* This function is called from run_a_cycle() every cycle */ 
  /* You do not add new code here */ 
  /* insert D-cache/I-cache (D-cache for only Lab #2) and wakes up instructions */ 
  fill_queue(); 

  /* move memory requests from dram to cache and MSHR*/ /* out queue */ 
  send_bus_out_queue(); 

  /* memory requests are scheduled */ 
  dram_schedule(); 

  /* memory requests are moved from bus_queue to DRAM scheduler */
  push_dram_sch_queue();
  
  /* new memory requests send from MSRH to in_bus_queue */ 
  send_bus_in_queue(); 

}

/*******************************************************************/
/* get a new mshr entry 
/*******************************************************************/

m_mshr_entry_s* memory_c::allocate_new_mshr_entry(void)
{
  
  if (m_mshr_free_list.empty())
    return NULL; 
  
  m_mshr_entry_s* entry = m_mshr_free_list.back(); 
  m_mshr_free_list.pop_back(); 
  m_mshr.push_back(entry); 
  
  return entry; 
  
 }






/* insert memory request into the MSHR, if there is no space return false */ 
/*******************************************************************/
/* memory_c::insert_mshr
/*******************************************************************/
bool memory_c::insert_mshr(Op *mem_op)
{

  bool insert=false; 
 
  //  step 1. create a new memory request 
  m_mshr_entry_s * entry = allocate_new_mshr_entry();
 
  //  step 2. no free entry means no space return ; 
  if (!entry) return false; // cannot insert a request into the mshr; 

  // insert op to into the mshr 
  entry->req_ops.push_back(mem_op);

    // step 3. initialize the memory request here 
  mem_req_s *t_mem_req = entry->m_mem_req; 
  
  if (mem_op->mem_type == MEM_LD ) { 
    t_mem_req->m_type = MRT_DFETCH; 
    t_mem_req->m_addr = mem_op->ld_vaddr; 
    t_mem_req->m_size = mem_op->mem_read_size;
  }
  else {
    t_mem_req->m_type = MRT_DSTORE; 
    t_mem_req->m_addr = mem_op->st_vaddr; 
    t_mem_req->m_size = mem_op->mem_write_size;
  }
    t_mem_req->m_rdy_cycle = 0; 
    t_mem_req->m_dirty = false; 
    t_mem_req->m_done = false; 
    t_mem_req->m_id = ++m_unique_m_count; 
    t_mem_req->m_state = MEM_NEW; 
    
    entry->valid  = true; 
    entry->insert_time = cycle_count; 
    
    return true; 
}


/*******************************************************************/
/* searching for matching mshr entry using memory address 
/*******************************************************************/


m_mshr_entry_s * memory_c::search_matching_mshr(ADDRINT addr) {
  
  list<m_mshr_entry_s *>::const_iterator cii; 
  
  for (cii= m_mshr.begin() ; cii != m_mshr.end(); cii++) {
    
    m_mshr_entry_s* entry = (*cii);
    if (entry->valid) {
      mem_req_s *req = entry->m_mem_req; 
      if (req) {
	if ((req->m_addr)/m_block_size == (addr/m_block_size)) 
	  return entry; 
      }
    }
  }
  return NULL;
}

/*******************************************************************/
/* searching for matching mshr entry using memory address and return the iterator 
/*******************************************************************/

list<m_mshr_entry_s*>::iterator memory_c::search_matching_mshr_itr(ADDRINT addr) {
  
  list<m_mshr_entry_s *>::iterator cii; 
  
  for (cii= m_mshr.begin() ; cii != m_mshr.end(); cii++) {
    
    m_mshr_entry_s* entry = (*cii);
    if (entry->valid) {
      mem_req_s *req = entry->m_mem_req; 
      if (req) {
	if ((req->m_addr)/m_block_size == (addr/m_block_size)) 
	  return cii; 
      }
    }
  }
  return cii;
}



/*******************************************************************/
/*  search MSHR entries and find a matching MSHR entry and piggyback 
/*******************************************************************/

bool memory_c::check_piggyback(Op *mem_op)
{
  bool match = false; 
  ADDRINT addr;
  if (mem_op->mem_type == MEM_LD) 
    addr = mem_op->ld_vaddr;
  else 
    addr = mem_op->st_vaddr; 
  
  m_mshr_entry_s *entry = search_matching_mshr(addr);
  
  if(!entry)
  {
    return false;
  } else 
  {
    entry->req_ops.push_back(mem_op); 
    return true;
  }
}



/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*****    COMPLETE the following functions               ***********/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/




/*******************************************************************/
/* send a request from mshr into in_queue 
/*******************************************************************/

void memory_c::send_bus_in_queue() 
{
 
   // Traverse through the MSHR entries. If an entry in MSHR is not valid, continue through the loop. If valid, 
   // check if its a new entry in MSHR. If its a new entry, insert new memory request into dram_in_queue
	
	list <m_mshr_entry_s *>::const_iterator cii; 
	  
	for (cii= m_mshr.begin() ; cii != m_mshr.end(); cii++) 
	{
	    
		m_mshr_entry_s* entry = (*cii);
	    if (!entry->valid) 
	    	continue;
	    
	    mem_req_s* m_req=entry->m_mem_req;
	    
	    if(m_req->m_state==MEM_NEW) 
	    {  
	    	//insert mem request into dram_in_queue and update the memory state for the mem_req
	    	dram_in_queue.push_back(m_req);
	    	m_req->m_state=MEM_DRAM_IN;    
	    	return;
	    }	
	    
	}
	    	 
  // *END**** This is for the TA
} 



/*******************************************************************/
/* search from dram_schedule queue and scheule a request 
/*******************************************************************/

 void memory_c::dram_schedule() 
 {
   // Traverse all entries in all the dram_bank_sch_queues. 
   // If there is an unscheduled memory request, check the corresponding DRAM bank and see whether the DRAM bank is available. 
   // If the corresponding bank is idle, check the last row buffer number.
   // If the memory request has the same row id, its a row buffer hit. So, you set req->m_rdy_cycle = cycle_count + KNOB_MEM_LATENCY_ROW_HIT
   // Else if row buffer miss, you set req->m_rdy_cycle = cycle_count + KNOB_MEM_LATENCY_ROW_MISS, and also set row buffer id as the row id of the memory request. 
   // You also set dram_bank_rdy_cycle = req->m_rdy_cycle;
   
	 int req_dram_row_id;
	 
	 list<mem_req_s *>::const_iterator cii; 
	 
	 for(int i=0;i<m_dram_bank_num;i++)
	 {
		if(dram_bank_sch_queue[i].empty())
		 	continue;
		 
		for (cii= dram_bank_sch_queue[i].begin() ; cii != dram_bank_sch_queue[i].end(); cii++) 
		 {
		   mem_req_s* m_req = (*cii);   
		 
		   if(m_req->m_state==MEM_DRAM_SCH) // unscheduled memory request waiting to be scheduled
		    {   
			 if (dram_bank_rdy_cycle[i] < cycle_count) // DRAM bank output is available in row buffer
			 {
				req_dram_row_id=get_dram_row_id(m_req->m_addr); 
				
				if(dram_bank_open_row[i]==req_dram_row_id) // ASK: is this right way of comparing??? row buffer hit
				  {
					m_req->m_rdy_cycle = cycle_count + KNOB(KNOB_MEM_LATENCY_ROW_HIT)->getValue();
					dram_row_buffer_hit_count++;
				  }
				   
				else	
				   {
					dram_row_buffer_miss_count++;
					m_req->m_rdy_cycle = cycle_count + KNOB(KNOB_MEM_LATENCY_ROW_MISS)->getValue();
					
				   }
				   
				dram_bank_open_row[i]=req_dram_row_id;
				dram_bank_rdy_cycle[i] = m_req->m_rdy_cycle;
				m_req->m_state=MEM_DRAM_DONE;
				
				break; // since only 1 mem request in a particular bank can be scheduled
			 }
			 else
				 break; // dont iterate through the rest of the mem requests since cycle_count is < dram_bank_rdy_cycle[i] for all of them.
		    }
		} 
	 }
}


/*******************************************************************/
/*  push_dram_sch_queue()
/*******************************************************************/

void memory_c::push_dram_sch_queue()
{
  // If dram_in_queue is empty, return. Else get the first memory request in dram_in_queue,  
  // find out the dram_bank_id and move a memory request from the dram_in_queue to the corresponding 
  // dram_bank_sch_queue  
    
	int dram_bank_id;
	if(dram_in_queue.empty())
	{ 
		return;}
	
	mem_req_s* m_req=dram_in_queue.front();
	dram_bank_id=get_dram_bank_id(m_req->m_addr);
	
	dram_bank_sch_queue[dram_bank_id].push_back(m_req);
	m_req->m_state=MEM_DRAM_SCH;
	dram_in_queue.pop_front();
	
	   
}

/*******************************************************************/
/* send bus_out_queue 
/*******************************************************************/

void memory_c::send_bus_out_queue() 
{

   /* check dram_bank_sch_queue and make sure its not empty. 
    * When the memory request is ready, the memory request is moved into the dram_out_queue. 
      If multiple requests are ready at the same cycle, we still move all the requests into the dram_out_queue  */
	
	 for(int i=0;i<m_dram_bank_num;i++)
		 {
			 if(dram_bank_sch_queue[i].empty())
			 	continue;
			 
			 mem_req_s* m_req=dram_bank_sch_queue[i].front();
			 if(m_req->m_state==MEM_DRAM_DONE)
			 {	 
				 if (dram_bank_rdy_cycle[i] < cycle_count) // DRAM bank output is available in row buffer
				  {
				     // memory request is ready
					 m_req->m_state=MEM_DRAM_OUT; 
					 dram_out_queue.push_back(m_req);
					 dram_bank_sch_queue[i].pop_front();
					 
				  }
			 }
		 }
}	 
  



/*******************************************************************/
/*  fill_queue 
/*******************************************************************/

void memory_c::fill_queue() 
{

   /* For Lab #2, you need to fill out this function */ 
  /* CAUTION!: This function is not completed. Please complete this function */ 
  
  if (dram_out_queue.empty()) return; 
  
    mem_req_s *req = dram_out_queue.front(); 
    dcache_insert(req->m_addr); // insert into dcache
    
    dram_out_queue.pop_front(); 
    
    /* search for matching mshr  entry */ 
    m_mshr_entry_s *entry = search_matching_mshr(req->m_addr); 
    
    
    while(entry->req_ops.size()) 
    {
      Op *w_op = entry->req_ops.front();
      broadcast_rdy_op(w_op); 
      entry->req_ops.pop_front(); 
      
    }
    
    /* The following code will free mshr entry */ 
    list<m_mshr_entry_s *>::iterator mii = search_matching_mshr_itr(req->m_addr); 
    m_mshr.erase(mii); 
  
  free_mshr_entry(entry); 
  
}





/*******************************************************************/
/*  store-load forwarding features, cache addresses, cache load/store types 
/*******************************************************************/


bool memory_c::store_load_forwarding(Op *mem_op)
{

  /* This function is called only if mem_op is a LOAD instruction. Search for a matching address in the MSHR. If found, check 
	if forwarding is possible. If yes, return TRUE else return FALSE */
	    bool can_forward=false;
	    ADDRINT current_m_addr,current_m_size;
	    
	    if(mem_op->mem_type==MEM_LD) // for store-load forwarding
	    	{
	    	current_m_addr=mem_op->ld_vaddr;
	    	current_m_size=(ADDRINT)mem_op->mem_read_size;
	    	}
	    
	    
		m_mshr_entry_s *entry = search_matching_mshr(current_m_addr);
		
		if(entry!=NULL)
		{	
			
		  list <Op *>::const_iterator cii;
		  
		  for(cii=entry->req_ops.begin(); cii != entry->req_ops.end(); cii++)
		  {
			  Op* op=(*cii);
			  if (op->mem_type == MEM_ST)
			   { 
			        if ((op->st_vaddr <= current_m_addr) && (op->st_vaddr + op->mem_write_size >= current_m_addr + current_m_size)) 
			        {
			        	// store-load forwarding possible
			        	can_forward=true;
			        	break;
			        }
			        	
			   } 
		  }
		}	  
		  			
    return can_forward;

}

/*******************************************************************/
/*  store-store forwarding features, cache addresses, cache load/store types 
/*******************************************************************/

bool memory_c::store_store_forwarding(Op *mem_op)
{

  /* This function is called only if mem_op is a LOAD instruction. Search for a matching address in the MSHR. If found, check 
	if forwarding is possible. If yes, return TRUE else return FALSE */
	    bool can_forward=false;
	    ADDRINT current_m_addr,current_m_size;
	    
	 	    
	    if(mem_op->mem_type==MEM_ST) // for store-store forwarding
	   	    {
	    	current_m_addr=mem_op->st_vaddr;
	    	current_m_size=(ADDRINT)mem_op->mem_write_size;
	   	    }
	    
		m_mshr_entry_s *entry = search_matching_mshr(current_m_addr);
		
		if(entry!=NULL)
		{	
			
		  list <Op *>::const_iterator cii;
		  
		  for(cii=entry->req_ops.begin(); cii != entry->req_ops.end(); cii++)
		  {
			  Op* op=(*cii);
			  if (op->mem_type == MEM_ST)
			   { 
			        if ((op->st_vaddr <= current_m_addr) && (op->st_vaddr + op->mem_write_size >= current_m_addr + current_m_size)) 
			        {
			        	// store-store forwarding possible
			        	can_forward=true;
			        	break;
			        }
			        	
			   } 
		  }
		}	  
		  			
    return can_forward;

}

