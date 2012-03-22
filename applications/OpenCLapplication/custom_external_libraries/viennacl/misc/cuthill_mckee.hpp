#ifndef VIENNACL_MISC_CUTHILL_MCKEE_HPP
#define VIENNACL_MISC_CUTHILL_MCKEE_HPP

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


/** @file viennacl/misc/cuthill_mckee.hpp
*    @brief Implementation of several flavors of the Cuthill-McKee algorithm.  Experimental in 1.2.x.
*    
*   Contributed by Philipp Grabenweger, interface adjustments by Karl Rupp.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <deque>
#include <cmath>


namespace viennacl
{
  
  namespace detail
  {
    
    // function to calculate the increment of combination comb.
    // parameters:
    // comb: pointer to vector<int> of size m, m <= n
    //       1 <= comb[i] <= n for 0 <= i < m
    //       comb[i] < comb[i+1] for 0 <= i < m - 1
    //       comb represents an unordered selection of m values out of n
    // n: int
    //    total number of values out of which comb is taken as selection
    inline bool comb_inc(std::vector<int> & comb, int n)
    {
        int m;
        int k;
        
        m = comb.size();
        // calculate k as highest possible index such that (*comb)[k-1] can be incremented
        k = m;
        while ( (k > 0) && ( ((k == m) && (comb[k-1] == n)) || 
                            ((k < m) && (comb[k-1] == comb[k] - 1) )) )
        {
            k--;
        }
        if (k == 0) // no further increment of comb possible -> return false
        {
            return false;
        }
        else
        {
            (comb[k-1])++; // increment (*comb)[k-1],
            for (int i = k; i < m; i++) // and all higher index positions of comb are 
            // calculated just as directly following integer values (lowest possible values)
            {
                comb[i] = comb[k-1] + (i - k + 1);
            }
            return true;
        }
    }


    // function to generate a node layering as a tree structure rooted at
    // node s
    template <typename MatrixType>
    void generate_layering(MatrixType const & matrix, 
                           std::vector< std::vector<int> > & l,
                           int s)
    {
      std::size_t n = matrix.size();
      //std::vector< std::vector<int> > l;
      std::vector<bool> inr(n, false);
      std::vector<int> nlist;
      
      nlist.push_back(s);
      inr[s] = true;
      l.push_back(nlist);
      
      for (;;)
      {
          nlist.clear();
          for (std::vector<int>::iterator it  = l.back().begin(); 
                                          it != l.back().end();
                                          it++)
          {
              for (typename MatrixType::value_type::const_iterator it2  = matrix[*it].begin(); 
                                                         it2 != matrix[*it].end();
                                                         it2++)
              {
                  if (it2->first == *it) continue;
                  if (inr[it2->first]) continue;
                  
                  nlist.push_back(it2->first);
                  inr[it2->first] = true;
              }
          }
          
          if (nlist.size() == 0)
              break;

          l.push_back(nlist);
      }
      
    }

    
    // comparison function for comparing two vector<int> values by their 
    // [1]-element
    inline bool cuthill_mckee_comp_func(std::vector<int> const & a,
                                        std::vector<int> const & b)
    {
        return (a[1] < b[1]);
    }
    
  }
  
  //
  // Part 1: The original Cuthill-McKee algorithm
  //
  
  
  struct cuthill_mckee_tag {};
  
  /** @brief Function for the calculation of a node number permutation to reduce the bandwidth of an incidence matrix by the Cuthill-McKee algorithm
   * 
   * references:
   *    Algorithm was implemented similary as described in 
   *      "Tutorial: Bandwidth Reduction - The CutHill-
   *      McKee Algorithm" posted by Ciprian Zavoianu as weblog at
   *    http://ciprian-zavoianu.blogspot.com/2009/01/project-bandwidth-reduction.html
   *    on January 15, 2009
   *    (URL taken on June 14, 2011) 
   * 
   * @param matrix  vector of n matrix rows, where each row is a map<int, double> containing only the nonzero elements
   * @return permutation vector r. r[l] = i means that the new label of node i will be l.
   *
   */
  template <typename MatrixType>
  std::vector<int> reorder(MatrixType const & matrix, cuthill_mckee_tag)
  {
    std::size_t n = matrix.size();
    std::vector<int> r;
    std::vector<bool> inr(n, false); // status array which remembers which nodes have been added to r
    std::deque<int> q;
    std::vector< std::vector<int> > nodes;
    std::vector<int> tmp(2);
    int p = 0;
    int c;
    
    int deg;
    int deg_min;
    
    r.reserve(n);
    nodes.reserve(n);
    
    do
    {
        // under all nodes not yet in r determine one with minimal degree
        deg_min = -1;
        for (std::size_t i = 0; i < n; i++)
        {
            if (!inr[i])
            {
                deg = matrix[i].size() - 1; // node degree
                if (deg_min < 0 || deg < deg_min)
                {
                    p = i; // node number
                    deg_min = deg;
                }
            }
        }
        q.push_back(p); // push that node p with minimal degree on q
        
        do
        {
            c = q.front();
            q.pop_front();
            if (!inr[c])
            {
                r.push_back(c);
                inr[c] = true;
                
                // add all neighbouring nodes of c which are not yet in r to nodes
                nodes.resize(0);
                for (typename MatrixType::value_type::const_iterator it =  matrix[c].begin(); it != matrix[c].end(); it++)
                {
                    if (it->first == c) continue;
                    if (inr[it->first]) continue;
                    
                    tmp[0] = it->first;
                    tmp[1] = matrix[it->first].size() - 1;
                    nodes.push_back(tmp);
                }
                
                // sort nodes by node degree
                std::sort(nodes.begin(), nodes.end(), detail::cuthill_mckee_comp_func);
                for (std::vector< std::vector<int> >::iterator it = nodes.begin(); it != nodes.end(); it++)
                {
                    q.push_back((*it)[0]);
                }
            }
        } while (q.size() != 0);
    } while (r.size() < n);
    
    return r;
  }
  
  
  //
  // Part 2: Advanced Cuthill McKee
  //

  /** @brief Tag for the advanced Cuthill-McKee algorithm */ 
  class advanced_cuthill_mckee_tag
  {
    public:
      /** @brief CTOR which may take the additional parameters for the advanced algorithm.
        * 
        * additional parameters for CTOR:
        *   a:  0 <= a <= 1
        *     parameter which specifies which nodes are tried as starting nodes
        *     of generated node layering (tree structure whith one ore more 
        *     starting nodes).
        *     the relation deg_min <= deg <= deg_min + a * (deg_max - deg_min) 
        *     must hold for node degree deg for a starting node, where deg_min/
        *     deg_max is the minimal/maximal node degree of all yet unnumbered
        *     nodes.
        *    gmax:
        *      integer which specifies maximum number of nodes in the root
        *      layer of the tree structure (gmax = 0 means no limit)
        * 
        * @return permutation vector r. r[l] = i means that the new label of node i will be l.
        *
       */
      advanced_cuthill_mckee_tag(double a = 0.0, std::size_t gmax = 1) : starting_node_param_(a), max_root_nodes_(gmax) {}
      
      double starting_node_param() const { return starting_node_param_;}
      void starting_node_param(double a) { if (a >= 0) starting_node_param_ = a; }
      
      std::size_t max_root_nodes() const { return max_root_nodes_; }
      void max_root_nodes(std::size_t gmax) { max_root_nodes_ = gmax; }      
      
    private:
      double starting_node_param_;
      std::size_t max_root_nodes_;
  };
  


  /** @brief Function for the calculation of a node number permutation to reduce the bandwidth of an incidence matrix by the advanced Cuthill-McKee algorithm
   * 
   *
   *  references:
   *    see description of original Cuthill McKee implementation, and
   *    E. Cuthill and J. McKee: "Reducing the Bandwidth of sparse symmetric Matrices".
   *    Naval Ship Research and Development Center, Washington, D. C., 20007
   */
  template <typename MatrixType>
  std::vector<int> reorder(MatrixType const & matrix,
                           advanced_cuthill_mckee_tag const & tag)
  {
    std::size_t n = matrix.size();
    double a = tag.starting_node_param();
    std::size_t gmax = tag.max_root_nodes();
    std::vector<int> r;
    std::vector<int> r_tmp;
    std::vector<int> r_best;
    std::vector<int> r2(n);
    std::vector<bool> inr(n, false);
    std::vector<bool> inr_tmp(n);
    std::vector<bool> inr_best(n);
    std::deque<int> q;
    std::vector< std::vector<int> > nodes;
    std::vector<int> nodes_p;
    std::vector<int> tmp(2);
    std::vector< std::vector<int> > l;
    int deg_min;
    int deg_max;
    int deg_a;
    int deg;
    int bw;
    int bw_best;
    std::vector<int> comb;
    std::size_t g;
    int c;
    
    r.reserve(n);
    r_tmp.reserve(n);
    r_best.reserve(n);
    nodes.reserve(n);
    nodes_p.reserve(n);
    comb.reserve(n);
    
    do
    {   
        // add to nodes_p all nodes not yet in r which are candidates for the root node layer  
        // search unnumbered node and generate layering 
        for (std::size_t i = 0; i < n; i++)
        {
            if (!inr[i])
            {
                detail::generate_layering(matrix, l, i);
                break;
            }
        }
        nodes.resize(0);
        for (std::vector< std::vector<int> >::iterator it = l.begin();
          it != l.end(); it++)
        {
            for (std::vector<int>::iterator it2 = it->begin();
              it2 != it->end(); it2++)
            {
                tmp[0] = *it2;
                tmp[1] = matrix[*it2].size() - 1;
                nodes.push_back(tmp);
            }
        }
        // determine minimum and maximum node degree
        deg_min = -1;
        deg_max = -1;
        for (std::vector< std::vector<int> >::iterator it = nodes.begin(); 
          it != nodes.end(); it++)
        {
            deg = (*it)[1];
            if (deg_min < 0 || deg < deg_min)
            {
                deg_min = deg;
            }
            if (deg_max < 0 || deg > deg_max)
            {
                deg_max = deg;
            }
        }
        deg_a = deg_min + (int) (a * (deg_max - deg_min));
        nodes_p.resize(0);
        for (std::vector< std::vector<int> >::iterator it = nodes.begin(); 
          it != nodes.end(); it++)
        {
            if ((*it)[1] <= deg_a)
            {
                nodes_p.push_back((*it)[0]);
            }
        }
        
        inr_tmp = inr;
        g = 1;
        comb.resize(1);
        comb[0] = 1;
        bw_best = -1;
        
        for (;;) // for all combinations of g <= gmax root nodes repeat
        {
            inr = inr_tmp;
            r_tmp.resize(0);
            
            // add the selected root nodes according to actual combination comb to q
            for (std::vector<int>::iterator it = comb.begin(); 
              it != comb.end(); it++)
            {
                q.push_back(nodes_p[(*it)-1]);
            }
  
            do // perform normal CutHill-McKee algorithm for given root nodes with 
            // resulting numbering stored in r_tmp
            {
                c = q.front();
                q.pop_front();
                if (!inr[c])
                {
                    r_tmp.push_back(c);
                    inr[c] = true;
                    
                    nodes.resize(0);
                    for (typename MatrixType::value_type::const_iterator it = matrix[c].begin(); it != matrix[c].end(); it++)
                    {
                        if (it->first == c) continue;
                        if (inr[it->first]) continue;
                        
                        tmp[0] = it->first;
                        tmp[1] = matrix[it->first].size() - 1;
                        nodes.push_back(tmp);
                    }
                    std::sort(nodes.begin(), nodes.end(), detail::cuthill_mckee_comp_func);
                    for (std::vector< std::vector<int> >::iterator it = 
                      nodes.begin(); it != nodes.end(); it++)
                    {
                        q.push_back((*it)[0]);
                    }
                }
            } while (q.size() != 0);
            
            // calculate resulting bandwith for root node combination
            // comb for current numbered component of the node graph
            for (std::size_t i = 0; i < r_tmp.size(); i++)
            {
                r2[r_tmp[i]] = r.size() + i;
            }
            bw = 0;
            for (std::size_t i = 0; i < r_tmp.size(); i++)
            {
                for (typename MatrixType::value_type::const_iterator it  = matrix[r_tmp[i]].begin(); 
                                                                     it != matrix[r_tmp[i]].end();
                                                                     it++)
                {
                    bw = std::max(bw, std::abs(static_cast<int>(r.size() + i) - r2[it->first]));
                }
            }
            
            // remember ordering r_tmp in r_best for smallest bandwith
            if (bw_best < 0 || bw < bw_best)
            {
                r_best = r_tmp;
                bw_best = bw;
                inr_best = inr;
            }
            
            // calculate next combination comb, if not existing
            // increment g if g stays <= gmax, or else terminate loop
            if (!detail::comb_inc(comb, nodes_p.size()))
            {
                g++;
                if ( (gmax > 0 && g > gmax) || g > nodes_p.size())
                {
                    break;
                }
                comb.resize(g);
                for (std::size_t i = 0; i < g; i++)
                {
                    comb[i] = i + 1;
                }
            }
        }
        
        // store best order r_best in result array r
        for (std::vector<int>::iterator it = r_best.begin(); 
          it != r_best.end(); it++)
        {
            r.push_back((*it));
        }
        inr = inr_best;
        
    } while (r.size() < n);
    
    return r;
  }
  
} //namespace viennacl
    

#endif
