


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


   

   /***********************************************************************************************
    * LinkedListNode()
    * 
    * Default Constructor
    * 
    * Parameters:
    * Returns:
    ***********************************************************************************************/
    LinkedListNode::LinkedListNode() : data_(nullptr), next_(nullptr) {} 




   /***********************************************************************************************
    * LinkedListNode(T*)
    * 
    * Constructs LinkedListNode and sets data pointer
    * 
    * Template Parameters:
    *     T             Data type
    * Parameters:
    *     data          Data
    * Returns:
    ***********************************************************************************************/
    template<class T>
    LinkedListNode::LinkedListNode(T* data) : data_((void*)data), next_(nullptr) {}




   /***********************************************************************************************
    * void set_data(T*)
    * 
    * Sets pointer to data
    * 
    * Template Parameters:
    *     T             Data type
    * Parameters:
    *     data          Data
    * Returns:
    ***********************************************************************************************/ 
    template<class T>
    void LinkedListNode::set_data(T* data)
    {this->data_ = (void*)data;}


   /***********************************************************************************************
    * void* get_data()
    * 
    * Gets pointer to data
    * 
    * Parameters:
    * Returns:
    *     void*      pointer to data
    ***********************************************************************************************/ 
    void* LinkedListNode::get_data()
    {return this->data_;} 



   /***********************************************************************************************
    * void set_next_node(LinkedListNode* next)
    * 
    * Sets pointer to next node
    * 
    * Parameters:
    *     next          Next node 
    * Returns:
    ***********************************************************************************************/ 
    void LinkedListNode::set_next_node(LinkedListNode* next)
    {this->data_ = next;}



   /***********************************************************************************************
    * LinkedListNode* get_next_node()
    * 
    * Gets pointer to next node
    * 
    * Parameters:
    *     next          Next node 
    * Returns:
    ***********************************************************************************************/ 
    LinkedListNode* LinkedListNode::get_next_node()
    {return this->next_;}
 
}    















