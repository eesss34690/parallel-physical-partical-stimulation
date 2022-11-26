/*
 *  Octree_omp.cpp
 *  SpatialTest Project
 *
 *  Octree_omp implementation.
 *
 *  Created by radix on 05/04/08.
 *  Copyright Mykola Konyk, <mykola@konyk.org>, 2008.
 *  
 *  This code is under Microsoft Reciprocal License (Ms-RL)
 *  Please see http://www.opensource.org/licenses/ms-rl.html
 *
 *  Important points about the license (from Ms-RL):
 *
 *  [A] For any file you distribute that contains code from the software (in source code or binary format), you must provide 
 *  recipients the source code to that file along with a copy of this license, which license will govern that file. 
 *  You may license other files that are entirely your own work and do not contain code from the software under any terms 
 *  you choose.
 *
 *  [B] No Trademark License- This license does not grant you rights to use any contributors' name, logo, or trademarks.
 *
 *  [C] If you bring a patent claim against any contributor over patents that you claim are infringed by the software, your 
 *  patent license from such contributor to the software ends automatically.
 *
 *  [D] If you distribute any portion of the software, you must retain all copyright, patent, trademark, and attribution notices 
 *  that are present in the software.
 *
 *  [E] If you distribute any portion of the software in source code form, you may do so only under this license by including a
 *  complete copy of this license with your distribution. If you distribute any portion of the software in compiled or object 
 *  code form, you may only do so under a license that complies with this license.
 *
 *  [F] The software is licensed "as-is." You bear the risk of using it. The contributors give no express warranties, guarantees 
 *  or conditions. You may have additional consumer rights under your local laws which this license cannot change. To the extent 
 *  permitted under your local laws, the contributors exclude the implied warranties of merchantability, fitness for a particular 
 *  purpose and non-infringement.
 *
 */

#include "Base.h"
#include "Octree_omp.h"
#include <cstdlib>
#include <omp.h>

namespace SpatialTest
{
    //--
    int Octree_omp::s_i32MaxDepth                   = 5;
    float Octree_omp::s_f32MaxObjectNodeRatio       = 8.0f;

    //--
    Octree_omp::Octree_omp(const Vector3& refCenter, float f32HalfWidth, int i32Rebuild) :
        ISpatialStructure(),
        m_i32Rebuild(i32Rebuild),
        m_f32HalfWidth(f32HalfWidth)
    {
        // [rad] Create root node (0 depth)
        m_pRootNode = new Octree_ompNode(NULL, refCenter, f32HalfWidth, 0);
        
    }

    //--
    Octree_omp::~Octree_omp()
    {
        if(m_pRootNode)
        {
            delete(m_pRootNode);
            m_pRootNode = NULL;
        }
    }

    //--
    inline
    void
    Octree_omp::Preallocate(int i32Depth)
    {
        m_pRootNode->Preallocate(i32Depth);
    }

    //--
    void
    Octree_omp::VUpdate()
    {   
        ISpatialObject* pElement;
        
        // [rad] We are rebuilding octree every frame
        if(m_i32Rebuild)
        {
            ISpatialObject* pObjectList = NULL;
            
            // [rad] We are going to rebuild octree (without dealocating space)
            m_pRootNode->Rebuild();
        
        
            // [rad] Populate the octree
            std::vector<ISpatialObject*>::iterator iter_object;
            #pragma parallel omp for
            for(iter_object = m_vecObjects.begin(); iter_object != m_vecObjects.end(); iter_object++)
            {
                pElement = (*iter_object);
                
                pElement->VSetNext(pObjectList);
                pObjectList = pElement;
            }
            
            m_pRootNode->AddObjects(pObjectList, m_vecObjects.size());
        }
        else
        {      
        
            // [rad] Instead of rebuilding octree every frame we will remove objects
            // bottom-up
            Octree_ompNode* pNode;
            ISpatialObject* pObject;
            ISpatialObject* pPrev;
        
            std::vector<ISpatialObject*>::iterator iter_object;
            #pragma parallel omp for
            for(iter_object = m_vecObjects.begin(); iter_object != m_vecObjects.end(); iter_object++)
            {
                pElement = (*iter_object);
                pNode = static_cast<Octree_ompNode*>(pElement->VGetCell());
            
                // [rad] Remove element
                if(pNode->m_pObjects == pElement)
                {
                    // [rad] Remove object from the list (first position)
                    pNode->m_pObjects = pNode->m_pObjects->VGetNext();
                }
                else
                {
                    // [rad] traverse list and remove
                    pObject = pNode->m_pObjects;
                    #pragma omp parallel
                    while(pObject)
                    {
                        pPrev = pObject;
                        pObject = pObject->VGetNext();
                        
                        if(pObject == pElement)
                        {
                            // [rad] Remove object from the list
                            pPrev->VSetNext(pObject->VGetNext());
                            
                            break;
                        }
                    }
                }
                


                // [rad] Decrement node's object count
                pNode->m_i32ObjectCount--;
            
            
    
                // [rad] Attempt to re-insert element
                if(pNode->CheckContains(pElement))
                {
                    pNode->AddObject(pElement);
                }
                else
                {
                    // [rad] Go up one node
                    pNode = pNode->m_pParent;
                
                    // [rad] Check if this node still contains this object
                    #pragma omp parallel
                    while(pNode)
                    {
                        // [rad] Check if contains
                        if(pNode->CheckContains(pElement))
                        {
                            pNode->AddObject(pElement);
                            break;
                        }
                        else
                        {
                            pNode = pNode->m_pParent;
                        }
                    }
                }
            }
        }
        
        
        
        // [rad] Check collisionsgdb1
        std::vector<Octree_ompNode*> vecAncestors;
        m_pRootNode->CheckCollisions(vecAncestors);
        
    }

    //--
    void
    Octree_omp::VAddObjects(const std::vector<ISpatialObject*>& refObjects)
    {        
        float f32MinDiameter = 1.0e38f;
        float f32CellSize = 2.0f * m_f32HalfWidth;
        float f32Diameter;

        int i32Depth = 0;
        int i32Divisions = 1;

        std::vector<ISpatialObject*>::const_iterator iter_object;
        #pragma omp parallel for
        for(iter_object = refObjects.begin(); iter_object != refObjects.end(); iter_object++)
        {
            f32Diameter = (*iter_object)->VGetRadius();
            
            if(f32Diameter < f32MinDiameter)
            {
                // [rad] Found a smaller object, update
                f32MinDiameter = f32Diameter;
            }
        }
        
        // [rad] Make cell size 'x' times as big as the smallest diameter 
        f32MinDiameter *= s_f32MaxObjectNodeRatio;
        
        
        // [rad] Calculate depth
        f32CellSize = 2.0f * m_f32HalfWidth;
        while(f32MinDiameter <= f32CellSize)
        {
            // [rad] Don't make more than max depth 
            if(s_i32MaxDepth == i32Depth)
            {
                break;
            }
            
            i32Divisions *= 2;                        
            f32CellSize = 2.0f * m_f32HalfWidth / static_cast<float>(i32Divisions);        
            
            i32Depth++;
        }
        
        
        // [rad] Delete everything but the root
        m_pRootNode->Free();
        
        // [rad] Pre-allocate the tree
        Preallocate(i32Depth);


        ISpatialObject* pObjectList = NULL;
        ISpatialObject* pObjectTemp;

        
        // [rad] Insert the elements into the tree
        #pragma omp parallel for
        for(iter_object = refObjects.begin(); iter_object != refObjects.end(); iter_object++)
        {
            if(m_i32Rebuild)
            {
                m_pRootNode->AddObject((*iter_object));
            }
            else
            {
                pObjectTemp = (*iter_object);
            
                pObjectTemp->VSetNext(pObjectList);
                pObjectList = pObjectTemp;
            }
            
            
            // [rad] Store locally in a list, for fast iteration
            #pragma omp critical
            m_vecObjects.push_back((*iter_object));
        }
        
        
        if(!m_i32Rebuild)
        {
            // [rad] Add objects
            m_pRootNode->AddObjects(pObjectList, m_vecObjects.size());
        }
    }
    


    //--
    int Octree_ompNode::s_i32MinSplitCount          = 32;

    //--
    // **Octree_ompNode** constructor:
    // Construct the Octree_omp with all the child node nullptr and set other attributes
    Octree_ompNode::Octree_ompNode(Octree_ompNode* pParent, const Vector3& refCenter, float f32HalfWidth, int i32Depth) :
        m_vec3Center(refCenter),
        m_f32HalfWidth(f32HalfWidth),
        m_i32Depth(i32Depth),
        m_pObjects(NULL),
        m_pParent(pParent),
        m_i32ObjectCount(0)
    {
        for(int i = 0; i < 8; i++)
        {
            m_pChildren[i] = NULL;
        }
    }

    //--
    // **Octree_ompNode** Destructor:
    // Destruct the node and all of its children
    Octree_ompNode::~Octree_ompNode()
    {
        for(int i = 0; i < 8; i++)
        {
            if(m_pChildren[i])
            {
                delete(m_pChildren[i]);
                m_pChildren[i] = NULL;
            }
        }            
    }

    //--
    // **Octree_ompNode** Free:
    // [rad] Delete all children
    void
    Octree_ompNode::Free()
    {
        for(int i = 0; i < 8; i++)
        {
            if(m_pChildren[i])
            {
                delete(m_pChildren[i]);
                m_pChildren[i] = NULL;
            }
        }            
    }

    //--
    // **Octree_ompNode** CheckCollisions:
    // Iterate through all the objects
    // Check mutual collision first then iterate to the children
    void
    Octree_ompNode::CheckCollisions(std::vector<Octree_ompNode*>& refAncestors)
    {
        // [rad] Push this node
        refAncestors.push_back(this);
        
        // [rad] Check collisions between all objects
        std::vector<Octree_ompNode*>::iterator iter_node;
        for(iter_node = refAncestors.begin(); iter_node != refAncestors.end(); iter_node++)
        {
            (*iter_node)->CheckMutualCollisions(m_pObjects);
        }    
    
        // [rad] Recursively visit children
        for(int i = 0; i < 8; i++)
        {
            if(m_pChildren[i])
            {
                m_pChildren[i]->CheckCollisions(refAncestors);
            }
        }      
        
        // [rad] Pop this node
        refAncestors.pop_back();
    }

    //--
    // **Octree_ompNode** CheckContains:
    // Check the vicinity of the object and obtain the objects in range
    int
    Octree_ompNode::CheckContains(ISpatialObject* pObject)
    {
        int i32Index = 0;
        
        float f32Radius = pObject->VGetRadius();
        Vector3 vec3Center = pObject->VGetPosition(); 
        
        
        for(int i32Index = 0; i32Index < 3; i32Index++)
        {
            if((m_vec3Center[i32Index] - m_f32HalfWidth > vec3Center[i32Index] - f32Radius) ||
                (m_vec3Center[i32Index] + m_f32HalfWidth < vec3Center[i32Index] + f32Radius))
            {
                return(0);
            }
        }

        return(1);
    }

    //--
    // **Octree_ompNode** AddObjects:
    // Insert object into the list regarding to the number
    void
    Octree_ompNode::AddObjects(ISpatialObject* pObjectList, int i32ObjectCount)
    {
        ISpatialObject* pIter;
        ISpatialObject* pObject;
        
        // [rad] Check if we can stop splitting
        if(i32ObjectCount < s_i32MinSplitCount)
        {
            m_pObjects = pObjectList;
            
            pIter = m_pObjects;   
            while(pIter)
            {
                pIter->VSetCell(this);
                pIter = pIter->VGetNext();
            }
            
            m_i32ObjectCount = i32ObjectCount;
            
            return;
        }
        
        // add object start
        int i32Index = 0;
        ISpatialObject* apChildObjects[8];
        int aiChildCounts[8];
        
        for(i32Index = 0; i32Index < 8; i32Index++)
        {
            apChildObjects[i32Index] = NULL;
            aiChildCounts[i32Index] = 0;
        }
        
        
        pIter = pObjectList;
        #pragma omp parallel
        while(pIter)
        {
            pObject = pIter->VGetNext();
            
            int i32Position = 0;
            int i32Straddle = 0;
    
            float f32Radius = pIter->VGetRadius();
            Vector3 vec3Center = pIter->VGetPosition();
            
            for(i32Index = 0; i32Index < 3; i32Index++)
            {
                if(m_vec3Center[i32Index] < vec3Center[i32Index])
                {
                    if(m_vec3Center[i32Index] > vec3Center[i32Index] - f32Radius)
                    {
                        // [rad] Straddle occurs
                        i32Straddle = 1;
                        break;
                    }
                    else
                    {
                        i32Position |= (1 << i32Index);
                    }
                }
                else
                {
                    if(m_vec3Center[i32Index] < vec3Center[i32Index] + f32Radius)
                    {
                        // [rad] Straddle occurs
                        i32Straddle = 1;
                        break;
                    }
                }
            }
        
        
            if(!i32Straddle && m_pChildren[i32Position])
            {
                // [rad] Contained in existing child node
                //m_pChildren[i32Position]->AddObject(pObject);
                
                pIter->VSetNext(apChildObjects[i32Position]);
                apChildObjects[i32Position] = pIter;
                #pragma omp critical
                aiChildCounts[i32Position]++;
                
            }
            else
            {
                // [rad] Store this node for fast back-link
                pIter->VSetCell(this);
        
                // [rad] Straddling or no child node available
                pIter->VSetNext(m_pObjects);
                m_pObjects = pIter;
                #pragma omp critical
                m_i32ObjectCount++;
            }
            
            
            pIter = pObject;
        }
        
        
        // [rad] At this point
        #pragma omp parallel for
        for(i32Index = 0; i32Index < 8; i32Index++)
        {            
            // [rad] Delegate to children
            if(m_pChildren[i32Index] && aiChildCounts[i32Index])
            {
                m_pChildren[i32Index]->AddObjects(apChildObjects[i32Index], aiChildCounts[i32Index]);
            }
        }    
    }

    //--
    void
    Octree_ompNode::AddObject(ISpatialObject* pObject)
    {
        int i32Index = 0;
        int i32Position = 0;
        int i32Straddle = 0;
        
        float f32Radius = pObject->VGetRadius();
        Vector3 vec3Center = pObject->VGetPosition();      
        
        for(int i32Index = 0; i32Index < 3; i32Index++)
        {
            if(m_vec3Center[i32Index] < vec3Center[i32Index])
            {
                if(m_vec3Center[i32Index] > vec3Center[i32Index] - f32Radius)
                {
                    // [rad] Straddle occurs
                    i32Straddle = 1;
                    break;
                }
                else
                {
                    i32Position |= (1 << i32Index);
                }
            }
            else
            {
                if(m_vec3Center[i32Index] < vec3Center[i32Index] + f32Radius)
                {
                    // [rad] Straddle occurs
                    i32Straddle = 1;
                    break;
                }
            }
            
        }
        
        
        if(!i32Straddle && m_pChildren[i32Position])
        {
            // [rad] Contained in existing child node
            m_pChildren[i32Position]->AddObject(pObject);
        }
        else
        {
            // [rad] Store this node for fast back-link
            pObject->VSetCell(this);
        
            // [rad] Straddling or no child node available
            pObject->VSetNext(m_pObjects);
            m_pObjects = pObject;
            
            // [rad] Increment number of elements
            m_i32ObjectCount++;
        }

    }

    //--
    inline
    bool
    Octree_ompNode::IsLeaf() const
    {
        return(!m_pChildren[0] && !m_pChildren[1] && !m_pChildren[2] &&
                !m_pChildren[3] && !m_pChildren[4] && !m_pChildren[5] &&
                !m_pChildren[6] && !m_pChildren[7]);
    }

    //--
    void
    Octree_ompNode::Rebuild()
    {
        #pragma omp parallel for
        for(int i = 0; i < 8; i++)
        {
            if(m_pChildren[i])
            {
                m_pChildren[i]->Rebuild();
            }
            
            m_pObjects = NULL;
        }
    }
 
    //--
    void
    Octree_ompNode::Preallocate(int i32Depth)
    {
        // [rad] Check conditions
        if(i32Depth && IsLeaf())
        {
            Vector3 vec3NewCenter;
            
            float f32Step = m_f32HalfWidth * 0.5f;
            #pragma omp parallel for
            for(int i = 0; i < 8; i++)
            {
                vec3NewCenter = m_vec3Center;
                
                if(i & 1)
                {
                    vec3NewCenter.x += f32Step;
                }
                else
                {
                    vec3NewCenter.x -= f32Step;
                }
                
                
                if(i & 2)
                {
                    vec3NewCenter.y += f32Step;
                }
                else
                {
                    vec3NewCenter.y -= f32Step;
                }
                
                
                if(i & 4)
                {
                    vec3NewCenter.z += f32Step;
                }
                else
                {
                    vec3NewCenter.z -= f32Step;
                }
                
                
                // [rad] Allocate node
                m_pChildren[i] = new Octree_ompNode(this, vec3NewCenter, f32Step, m_i32Depth + 1);
                
                // [rad] Build node recursively
                m_pChildren[i]->Preallocate(i32Depth - 1);                
            }
        }
    }

    //--
    // **Octree_omp_ompNode** CheckMutualCollisions:
    // based on the same layer, check whether there is collision
    void
    Octree_ompNode::CheckMutualCollisions(ISpatialObject* pObject)
    {
        ISpatialObject* pIter1 = pObject;
        ISpatialObject* pIter2;
        
        while(pIter1)
        {
            pIter2 = m_pObjects;
            
            while(pIter2)
            {
                
                if(pIter1 == pIter2) 
                {
                    break;
                }
                
                
                if(pIter1->VCheckCollision(pIter2))
                {
                    // [rad] Mark both as in collision
                    pIter1->VCollisionOn();
                    pIter2->VCollisionOn();
                }
                
                pIter2 = pIter2->VGetNext();
            }
            
            pIter1 = pIter1->VGetNext();
        }
    }

}
