// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_
#include <bits/stdc++.h>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;
 private:
  struct node
  {
    Point<N> punto;
    ElemType valor; 
    node *left=0;
    node *right=0;
    node(const value_type &value)
    {
      punto = value.first;
      valor = value.second;
    }
    node(Point<N> pt, ElemType v)
    {
      punto = pt;
      valor = v;
    }
    node(Point<N> pt, ElemType v, node *r, node *l)
    {
      punto = pt;
      valor = v;
      right = r;
      left = l;
    }
  };
  node *root;
  size_t dimension_;
  size_t size_;
  bool find(node **&ptr, const Point<N> &pt)
  {
    int d = 0; //dimension 1
    for(ptr = &root; (*ptr)&&(*ptr)->punto!=pt;)
    {
      if(pt[d%dimension_] > (*ptr)->punto[d%dimension_])
        ptr = &((*ptr)->right);
      else
        ptr = &((*ptr)->left);
      d++;
    }
    return ((*ptr)!=0);
  }
  bool find(node **&ptr, const Point<N> &pt) const
  {
    int d = 0; //dimension 1
    for(ptr = &root; (*ptr)&&(*ptr)->punto!=pt;)
    {
      if(pt[d%dimension_] > (*ptr)->punto[d%dimension_])
        ptr = &((*ptr)->left);
      else
        ptr = &((*ptr)->right);
      d++;
    }
    return ((*ptr)!=0);
  }
  // node *copy_nodes(node *original)
  // {
  //   std::cout<<"entre"<<std::endl;
  //   if(!original){
  //     new_node = new node(original->punto, original->valor, copy
  //   }
  //   std::cout<<"entre"<<std::endl;
  //   return new_node;
  // }
  void knn_fun(Point<N> pt, node* r, int d, double radio, node *&cand) const;
};
template<size_t N, typename ElemType>
void KDTree<N, ElemType>::knn_fun(Point<N> pt, node *r, int d, double radio,
node *&cand) const
{
  if(!r)
    return;
  if(distance(r->punto, pt) <= radio)
  {
    radio = distance(r->punto, pt);
    cand = r;
  }
  int eje = d%dimension_;
  if(pt[eje] <= r->punto[eje])
  {
    knn_fun(pt, r->left, d++, radio, cand);
  }
  else{
    knn_fun(pt, r->left, d++, radio, cand);
  }
  if(fabs(r->punto[eje]-pt[eje]) <= cand->punto[eje])
  {
    knn_fun(pt, r->left, d++, radio, cand);
  }
  else
  {
    knn_fun(pt, r->right, d++, radio, cand);
  }
}
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  size_ = 0;
  dimension_ = N;
  root = nullptr;
}
template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  // TODO(me): Fill this in.
}
template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  KDTree tmp;
  root = (rhs.root, root);
  size_ = rhs.size_;
  dimension_ = rhs.dimension_;
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  return size_;
} 
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  if(size_ == 0)
    return true;
  return false;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
  node *ptr = root;
  int d = 0;
  for(;(ptr);)
  {
    if(ptr->punto == pt)
      return true;
    if(pt[d%dimension_] > ptr->punto[d%dimension_])
      ptr = ptr->right;
    else
      ptr = ptr->left;
    d++;
  }
  return false;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  node **ptr;
  bool j = find(ptr,pt);
  if(j)
    (*ptr)->valor = value;
  else
  {
    *ptr = new node(pt, value);
    size_++;
  }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
  node **ptr;
  if(!find(ptr,pt))
    insert(pt,0);
  ElemType& ob = (*ptr)->valor;
  return ob;
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
  node **ptr;
  bool j = find(ptr, pt);
  if(!j)
    throw std::out_of_range("");
  ElemType& ob = (*ptr)->valor;
  return ob;
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  bool j = false;
  node *ptr = root;
  int d = 0;
  for(;(ptr);)
  {
    if(ptr->punto == pt)
    {
      j = true;
      break;
    }
    if(pt[d%dimension_] > ptr->punto[d%dimension_])
      ptr = ptr->right;
    else
      ptr = ptr->left;
    d++;
  }
  if(!j)
    throw std::out_of_range("");
  ElemType& ob = ptr->valor;
  return ob;
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  node *current = root;
  node *cand;
  double radio = std::numeric_limits<double>::infinity();
  knn_fun(key,current, 0, radio, cand);
  // ElemType new_element;
  return cand->valor;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {

  std::vector<ElemType> values;
  return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
