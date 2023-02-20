


namespace matl
{
    class LinkedListNode
    {
    private:  
        void*           data_;
        LinkedListNode* next_;
    public:

        // Default constructor
        LinkedListNode();

        // Data Constructor
        template<class T>
        LinkedListNode(T* data);

        // Set pointer to data
        template<class T>
        void set_data(T* data);

        // Get pointer to data
        void* get_data();

        // Set address of next node
        void set_next_node(LinkedListNode* next);

        // Access next node
        LinkedListNode* get_next_node();
    };
 
}    















