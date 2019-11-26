#ifndef INDEX_MIN_PQ_H_
#define INDEX_MIN_PQ_H_

#include <algorithm>
#include <sstream>
#include <utility>
#include <vector>

template <typename K>
class IndexMinPQ {
 public:
  // Constructor with max number of indexes
  explicit IndexMinPQ(int capacity);
  // Return number of items
  unsigned int Size();
  // Return top (ie index associated to minimum key)
  unsigned int Top();
  // Remove top
  void Pop();
  // Associates @key with index @idx
  void Push(const K &key, unsigned int idx);
  // Return whether @idx is a valid index
  bool Contains(unsigned int idx);
  // Change key associated to index @idx
  void ChangeKey(const K &key, unsigned int idx);

 private:
  // Private members
  unsigned int capacity;
  unsigned int cur_size;
  std::vector<K> keys;
  std::vector<unsigned int> heap_to_idx;
  std::vector<unsigned int> idx_to_heap;

  // Helper methods for indices
  unsigned int Root() {
    return 1;
  }
  unsigned int Parent(unsigned int i) {
    return i / 2;
  }
  unsigned int LeftChild(unsigned int i) {
    return 2 * i;
  }
  unsigned int RightChild(unsigned int i) {
    return 2 * i + 1;
  }

  // Helper methods for node testing
  bool HasParent(unsigned int i) {
    return i != Root();
  }
  bool IsNode(unsigned int i) {
    return i <= cur_size;
  }
  bool GreaterNode(unsigned int i, unsigned int j) {
    // Return true if node at index i is greater than node at index j, false
    // otherwise
    return (keys[heap_to_idx[i]] > keys[heap_to_idx[j]]);
  }

  // Helper methods for restructuring
  void SwapNodes(unsigned int i, unsigned int j) {
    // Swap nodes in heap
    std::swap(heap_to_idx[i], heap_to_idx[j]);
    // Update inverse mappings
    idx_to_heap[heap_to_idx[i]] = i;
    idx_to_heap[heap_to_idx[j]] = j;
  }
  void PercolateUp(unsigned int i);
  void PercolateDown(unsigned int i);

  // Helper method to check heap-order (useful for debugging)
  void CheckHeapOrder(unsigned int i) {
    if (!IsNode(i))
      return;
    if (HasParent(i) && GreaterNode(Parent(i), i)) {
      std::stringstream ss;
      ss << "Heap order error: "
          << "Parent ("
            << Parent(i) << ": " << heap_to_idx[Parent(i)] << ", "
            << keys[heap_to_idx[Parent(i)]] << ")"
          << " bigger than Child ("
            << i << ": " << heap_to_idx[i] << ", "
            << keys[heap_to_idx[i]] << ")";
      throw std::runtime_error(ss.str());
    }
    CheckHeapOrder(LeftChild(i));
    CheckHeapOrder(RightChild(i));
  }
};

template <typename K>
IndexMinPQ<K>::IndexMinPQ(int capacity)
  : capacity(capacity),
    keys(capacity),
    heap_to_idx(capacity + 1),
    idx_to_heap(capacity, 0) {
      cur_size = 0;
    }

template <typename K>
unsigned int IndexMinPQ<K>::Size() {
  return cur_size;
}

template <typename K>
unsigned int IndexMinPQ<K>::Top(void) {
  if (!Size())
    throw std::underflow_error("Priority queue underflow!");

  return heap_to_idx[Root()];
}

template <typename K>
void IndexMinPQ<K>::PercolateUp(unsigned int i) {
  while (HasParent(i) && GreaterNode(Parent(i), i)) {
    SwapNodes(Parent(i), i);
    i = Parent(i);
  }
}

template <typename K>
void IndexMinPQ<K>::Push(const K &key, unsigned int idx) {
  if (idx >= capacity)
    throw std::overflow_error("Index invalid!");
  if (Contains(idx))
    throw std::runtime_error("Index already exists!");

  // todo: push key-value pair made of @key and @idx
  // 1. Insert item at the end
  //  - Set both mapping tables properly
  //  - Set key in key vector
  heap_to_idx[++cur_size] = idx;
  idx_to_heap[idx] = cur_size;
  keys[idx] = key;
  // 2. Percolate up
  PercolateUp(cur_size);
}

template <typename K>
void IndexMinPQ<K>::PercolateDown(unsigned int i) {
  // While node has at least one child (if one, necessarily on the left)
  while (IsNode(LeftChild(i))) {
    // Find smallest children between left and right if any
    unsigned int child = LeftChild(i);
    if (IsNode(RightChild(i)) && GreaterNode(LeftChild(i), RightChild(i)))
      child = RightChild(i);

    // Exchange node with child to restore heap-order if necessary
    if (GreaterNode(i, child))
      SwapNodes(i, child);
    else
      break;

    // Do it again, one level down
    i = child;
  }
}

template <typename K>
void IndexMinPQ<K>::Pop() {
  if (!Size())
    throw std::underflow_error("Empty priority queue!");

  // todo: remove min item
  // 3. Mark idx_to_heap
  idx_to_heap[heap_to_idx[Root()]] = 0;
  // 1. Move last item back to root and reduce heap's size
  heap_to_idx[Root()] = heap_to_idx[cur_size--];
  // 2. Restore heap order
  PercolateDown(Root());
  // (for debugging, check heap order)
}

template <typename K>
bool IndexMinPQ<K>::Contains(unsigned int idx) {
  if (idx >= capacity)
    throw std::overflow_error("Index invalid!");
  return (idx_to_heap[idx] != 0);
}

template <typename K>
void IndexMinPQ<K>::ChangeKey(const K &key, unsigned int idx) {
  if (idx >= capacity)
    throw std::overflow_error("Index invalid!");
  if (!Contains(idx))
    throw std::runtime_error("Index does not exist!");

  // todo: modify the key associated to index @idx
  // 1. Update key in key vector
  K old_key = keys[idx];
  keys[idx] = key;
  // 2. Restore heap-order
  if (key < old_key) {
    PercolateUp(idx_to_heap[idx]);
  } else if (key > old_key) {
    PercolateDown(idx_to_heap[idx]);
  } else {
    throw std::runtime_error("Key already inserted!");
  }
}

#endif  // INDEX_MIN_PQ_H_
