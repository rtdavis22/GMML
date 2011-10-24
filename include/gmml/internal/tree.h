#ifndef TREE_H
#define TREE_H

#include <stack>

namespace gmml {

template<typename T>
class Tree;

template<typename T>
class TreeNode {
  public:
    TreeNode *append_child(T data);

    TreeNode *parent() { return parent_; }
    TreeNode *first_child() { return first_child_; }
    TreeNode *next_sibling() { return next_sibling_; }

    operator T&() { return data_; }    

  private:
    TreeNode(T& data, TreeNode *parent) : data_(data), parent_(parent),
                                          first_child_(NULL),
                                          next_sibling_(NULL) {}
    //no copying or assignment allowed
    TreeNode(const TreeNode&);
    void operator=(const TreeNode&);

    T data_;
    TreeNode *parent_;
    TreeNode *first_child_;
    TreeNode *next_sibling_;

    friend class Tree<T>;
};

template<typename T>
class Tree {
  public:
    struct iterator {
        iterator() : ptr(NULL) {}
        iterator(TreeNode<T> *ptr) : ptr(ptr) {}

        //T *operator->() { return 
        T operator*() { return ptr->data_; }

        //TreeNode<T>& operator*() { return *ptr; }
        TreeNode<T>* operator->() { return ptr; }
        bool operator==(const iterator& rhs) { return ptr == rhs.ptr; }
        bool operator!=(const iterator& rhs) { return !(ptr == rhs.ptr); }

        TreeNode<T> *ptr;
    };

    Tree(T data) : root_(new TreeNode<T>(data, NULL)) {}
    ~Tree();

    Tree *clone();

    iterator root() { return iterator(root_); }

  private:
    Tree(TreeNode<T> *root) : root_(root) {}
    TreeNode<T> *root_;    
};

template<typename T>
inline Tree<T> *Tree<T>::clone() {
    Tree *new_tree = new Tree(*root_);
    if (!root_->first_child())
        return new_tree;
    std::stack<iterator> st;
    st.push(new_tree->root());
    iterator it = root();
    while (it->first_child()) {
        it = it->first_child();
        st.push(st.top()->append_child(*it));
    }
    st.pop();
    while (it != NULL) {
        if (it->next_sibling()) {
            it = it->next_sibling();
            st.push(st.top()->append_child(*it));
            while (it->first_child()) {
                it = it->first_child();
                st.push(st.top()->append_child(*it));
            }
            st.pop();
        }
        else {
            st.pop();
            it = it->parent();
        }
    }
    return new_tree;
}

template<typename T>
inline TreeNode<T> *TreeNode<T>::append_child(T new_data) {
    TreeNode<T> *new_node = new TreeNode<T>(new_data, this);
    if (first_child_ == NULL) {
        first_child_ = new_node;
        return new_node;
    }
    TreeNode<T> *it = first_child_;
    while (it->next_sibling_)
        it = it->next_sibling_;
    it->next_sibling_ = new_node;
    return new_node;
}

template<typename T>
struct PreOrderIterator : public Tree<T>::iterator {
    PreOrderIterator(typename Tree<T>::iterator it) 
            : Tree<T>::iterator(it) { 
        while (this->ptr->first_child()) 
            this->ptr = this->ptr->first_child(); 
    }
    PreOrderIterator& operator++() {
        if (this->ptr->next_sibling()) {
            this->ptr = this->ptr->next_sibling();
            while (this->ptr->first_child())
                this->ptr = this->ptr->first_child();
        }
        else if (this->ptr->parent()) {
            this->ptr = this->ptr->parent();
        }
        else {
            this->ptr = NULL;
        }
        return *this;
    }

    PreOrderIterator operator++(int) {
        PreOrderIterator it(*this);
        operator++();
        return it;
    }

    bool at_end() const { return this->ptr == NULL; }
};

template<typename T>
inline Tree<T>::~Tree() {
    PreOrderIterator<T> it = root();
    while (!it.at_end()) {
        TreeNode<T> *ptr = it.ptr;
        ++it;
        delete ptr;
    }
}

} //namespace gmml

#endif
