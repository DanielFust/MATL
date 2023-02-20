

namespace matl
{
  class HeteroGeneousLinkedList
  {
    private:
      LinkedListNode head_;
      index_t        size_;  
    public:
      HeteroGeneousLinkedList();

      ~HeteroGeneousLinkedList();

      // Add item to the end of the linked list
      template<class T>
      void add_item(T*);

      // Insert item at the specified location
      template<class T>
      void insert_item(T*, index_t);

      // Access data by index
      void* operator() (index_t);

      // Remove list item by index
      remove_item(index_t);

      // Remove list item by node address
      remove_item(LinkedListNode*);

      // Get number of items in list
      index_t size();
      
  };
}






