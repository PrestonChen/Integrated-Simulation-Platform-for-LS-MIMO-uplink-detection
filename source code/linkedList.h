/*
 * linkedList.h

 * Notice the first node of the linked list is the protect node
 * which can protect the first address of the linked list
 *  Created on: Dec 11, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef LINKEDLIST_H_
#define LINKEDLIST_H_
#include "gslCommon.h"
#define LEN sizeof(struct Node)
struct Node{
	int index;
	struct Node *next;
};

struct Node *create(int Nt){   //create linked list
	struct Node *head, *current,*tail;
	int count;
	for (count = 0; count < Nt; count++){
		  current = (struct Node*)calloc(1, LEN);
		  current->index = count;
		  current->next = NULL;
		if(count == 0){
          head = current;
          tail = current;
		}
		if(Nt == 1){
		  break;
		}
	    tail->next = current;
	    tail = current;
	}
	return (head);
}

struct Node *split(struct Node *index2_head, gsl_combination *subset){
	//split the complete linked list index1 into two sub linked list index1 and index2
	struct Node *index1_head;
	struct Node *current, *first, *Nulling, *index1_tail;
	current = index2_head;
	int count = 0;
	int N = gsl_combination_k(subset);
	while(current->next != NULL && count < N){
		if(current->next->index == gsl_combination_get(subset, count)){
			Nulling = current->next;
			current->next = (Nulling->next);
			Nulling->next = NULL;
			if(count == 0){
				index1_head = Nulling;
				index1_tail = Nulling;
			}else{
			index1_tail->next = Nulling;
			index1_tail = (index1_tail->next);
			}
			count++;
		}else{
		current = current->next;
		}
	}
	return (index1_head);
}


struct Node *splitArray(struct Node *index2_head, int *subset, int N){
	//split the complete linked list index1 into two sub linked list index1 and index2
	struct Node *index1_head;
	struct Node *current, *first, *Nulling, *index1_tail;
	current = index2_head;
	int count = 0;
	while(current->next != NULL && count < N){
		if(current->next->index == *(subset+count)){
			Nulling = current->next;
			current->next = (Nulling->next);
			Nulling->next = NULL;
			if(count == 0){
				index1_head = Nulling;
				index1_tail = Nulling;
			}else{
			index1_tail->next = Nulling;
			index1_tail = (index1_tail->next);
			}
			count++;
		}else{
		current = current->next;
		}
	}
	return (index1_head);
}

int get(int k, struct Node *head){   // get the index in the kth struct Node
	int index;
	int count = 0;
	struct Node *current = head;
	struct Node *Nulling = NULL;
	while(count < k){
		current = current->next;
		count++;
	}
	Nulling = current->next;
	index = Nulling->index;   //get the index
	current->next = Nulling->next;  //delete Nulling struct Node
	free(Nulling);
	Nulling = NULL;
	return index;
}


void freeLinkedList(struct Node *head){
	struct Node *currentNode;
	struct Node *nextNode;
	currentNode = head;
	while(currentNode != NULL){
		nextNode = currentNode->next;
		free(currentNode);
		currentNode = nextNode;
	}
}


#endif /* LINKEDLIST_H_ */
