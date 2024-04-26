/*
  * This is a port https://github.com/lohedges/aabbcc. The namespace has been changed from
  * 'aabb' to 'aabb_base'. The class name has been changed from 'AABB' to 'AABB_base'. Also
  * class 'Tree' has been renamed to 'Tree_base'. Changed std::vector<double> to std::array<double,3> for
  * better performance. Other modifactions are marked as such.
  * The original code included the following copyright notice:

  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2018 Lester Hedges <lester.hedges+aabbcc@gmail.com>

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.

  This code was adapted from parts of the Box2D Physics Engine,
  http://www.box2d.org
*/

#include "AABB_base.h"

namespace aabb_base
{
    AABB_base::AABB_base()
    {
    }

    AABB_base::AABB_base(unsigned int dimension)
    {
        assert(dimension >= 2);

        // lowerBound.resize(dimension);
        // upperBound.resize(dimension);
    }

    AABB_base::AABB_base(const std::array<double,3>& lowerBound_, const std::array<double,3>& upperBound_) :
        lowerBound(lowerBound_), upperBound(upperBound_)
    {
        // Validate the dimensionality of the bounds vectors.
        if (lowerBound.size() != upperBound.size())
        {
            throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        }

        // Validate that the upper bounds exceed the lower bounds.
        for (unsigned int i=0;i<lowerBound.size();i++)
        {
            // Validate the bound.
            if (lowerBound[i] > upperBound[i])
            {
                throw std::invalid_argument("[ERROR]: AABB_base lower bound is greater than the upper bound!");
            }
        }

        surfaceArea = computeSurfaceArea();
        centre = computeCentre();
    }

    double AABB_base::computeSurfaceArea() const
    {
        // Sum of "area" of all the sides.
        double sum = 0;

        // General formula for one side: hold one dimension constant
        // and multiply by all the other ones.
        for (unsigned int d1 = 0; d1 < lowerBound.size(); d1++)
        {
            // "Area" of current side.
            double product = 1;

            for (unsigned int d2 = 0; d2 < lowerBound.size(); d2++)
            {
                if (d1 == d2)
                    continue;

                double dx = upperBound[d2] - lowerBound[d2];
                product *= dx;
            }

            // Update the sum.
            sum += product;
        }

        return 2.0 * sum;
    }

    double AABB_base::getSurfaceArea() const
    {
        return surfaceArea;
    }

    void AABB_base::merge(const AABB_base& aabb1, const AABB_base& aabb2)
    {
        // assert(aabb1.lowerBound.size() == aabb2.lowerBound.size());
        // assert(aabb1.upperBound.size() == aabb2.upperBound.size());

        // lowerBound.resize(aabb1.lowerBound.size());
        // upperBound.resize(aabb1.lowerBound.size());

        for (unsigned int i=0;i<lowerBound.size();i++)
        {
            lowerBound[i] = std::min(aabb1.lowerBound[i], aabb2.lowerBound[i]);
            upperBound[i] = std::max(aabb1.upperBound[i], aabb2.upperBound[i]);
        }

        surfaceArea = computeSurfaceArea();
        centre = computeCentre();
    }

    bool AABB_base::contains(const AABB_base& aabb_base) const
    {
        assert(aabb_base.lowerBound.size() == lowerBound.size());

        for (unsigned int i=0;i<lowerBound.size();i++)
        {
            if (aabb_base.lowerBound[i] < lowerBound[i]) return false;
            if (aabb_base.upperBound[i] > upperBound[i]) return false;
        }

        return true;
    }

    bool AABB_base::overlaps(const AABB_base& aabb_base, bool touchIsOverlap) const
    {
        assert(aabb_base.lowerBound.size() == lowerBound.size());

        bool rv = true;

        if (touchIsOverlap)
        {
            for (unsigned int i = 0; i < lowerBound.size(); ++i)
            {
                if (aabb_base.upperBound[i] < lowerBound[i] || aabb_base.lowerBound[i] > upperBound[i])
                {
                    rv = false;
                    break;
                }
            }
        }
        else
        {
            for (unsigned int i = 0; i < lowerBound.size(); ++i)
            {
                if (aabb_base.upperBound[i] <= lowerBound[i] || aabb_base.lowerBound[i] >= upperBound[i])
                {
                    rv = false;
                    break;
                }
            }
        }

        return rv;
    }

    std::array<double,3> AABB_base::computeCentre()
    {
        std::array<double,3> position{};

        for (unsigned int i=0;i<position.size();i++)
            position[i] = 0.5 * (lowerBound[i] + upperBound[i]);

        return position;
    }

    void AABB_base::setDimension(unsigned int dimension)
    {
        assert(dimension >= 2);

        // lowerBound.resize(dimension);
        // upperBound.resize(dimension);
    }

    Node::Node()
    {
    }

    bool Node::isLeaf() const
    {
        return (left == NULL_NODE);
    }

    Tree_base::Tree_base(unsigned int dimension_,
               double skinThickness_,
               unsigned int nParticles,
               bool touchIsOverlap_) :
        dimension(dimension_), isPeriodic(false), skinThickness(skinThickness_),
        touchIsOverlap(touchIsOverlap_)
    {
        // Validate the dimensionality.
        if ((dimension < 2))
        {
            throw std::invalid_argument("[ERROR]: Invalid dimensionality!");
        }

        // Initialise the periodicity vector.
        periodicity.resize(dimension);
        std::fill(periodicity.begin(), periodicity.end(), false);

        // Initialise the tree.
        root = NULL_NODE;
        nodeCount = 0;
        nodeCapacity = nParticles;
        nodes.resize(nodeCapacity);

        // Build a linked list for the list of free nodes.
        for (unsigned int i=0;i<nodeCapacity-1;i++)
        {
            nodes[i].next = i + 1;
            nodes[i].height = -1;
        }
        nodes[nodeCapacity-1].next = NULL_NODE;
        nodes[nodeCapacity-1].height = -1;

        // Assign the index of the first free node.
        freeList = 0;
    }

    Tree_base::Tree_base(unsigned int dimension_,
               double skinThickness_,
               const std::vector<bool>& periodicity_,
               const std::array<double,3>& boxSize_,
               unsigned int nParticles,
               bool touchIsOverlap_) :
        dimension(dimension_), skinThickness(skinThickness_),
        periodicity(periodicity_), boxSize(boxSize_),
        touchIsOverlap(touchIsOverlap_)
    {
        // Validate the dimensionality.
        if (dimension < 2)
        {
            throw std::invalid_argument("[ERROR]: Invalid dimensionality!");
        }

        // Validate the dimensionality of the vectors.
        if ((periodicity.size() != dimension) || (boxSize.size() != dimension))
        {
            throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        }

        // Initialise the tree.
        root = NULL_NODE;
        nodeCount = 0;
        nodeCapacity = nParticles;
        nodes.resize(nodeCapacity);

        // Build a linked list for the list of free nodes.
        for (unsigned int i=0;i<nodeCapacity-1;i++)
        {
            nodes[i].next = i + 1;
            nodes[i].height = -1;
        }
        nodes[nodeCapacity-1].next = NULL_NODE;
        nodes[nodeCapacity-1].height = -1;

        // Assign the index of the first free node.
        freeList = 0;

        // Check periodicity.
        isPeriodic = false;
        // posMinImage.resize(dimension);
        // negMinImage.resize(dimension);
        for (unsigned int i=0;i<dimension;i++)
        {
            posMinImage[i] =  0.5*boxSize[i];
            negMinImage[i] = -0.5*boxSize[i];

            if (periodicity[i])
                isPeriodic = true;
        }
    }

    void Tree_base::setPeriodicity(const std::vector<bool>& periodicity_)
    {
        periodicity = periodicity_;
    }

    void Tree_base::setBoxSize(const std::array<double,3>& boxSize_)
    {
        boxSize = boxSize_;
    }

    unsigned int Tree_base::allocateNode()
    {
        // Exand the node pool as needed.
        if (freeList == NULL_NODE)
        {
            assert(nodeCount == nodeCapacity);

            // The free list is empty. Rebuild a bigger pool.
            nodeCapacity *= 2;
            nodes.resize(nodeCapacity);

            // Build a linked list for the list of free nodes.
            for (unsigned int i=nodeCount;i<nodeCapacity-1;i++)
            {
                nodes[i].next = i + 1;
                nodes[i].height = -1;
            }
            nodes[nodeCapacity-1].next = NULL_NODE;
            nodes[nodeCapacity-1].height = -1;

            // Assign the index of the first free node.
            freeList = nodeCount;
        }

        // Peel a node off the free list.
        unsigned int node = freeList;
        freeList = nodes[node].next;
        nodes[node].parent = NULL_NODE;
        nodes[node].left = NULL_NODE;
        nodes[node].right = NULL_NODE;
        nodes[node].height = 0;
        nodes[node].aabb_base.setDimension(dimension);
        nodeCount++;

        return node;
    }

    void Tree_base::freeNode(unsigned int node)
    {
        assert(node < nodeCapacity);
        assert(0 < nodeCount);

        nodes[node].next = freeList;
        nodes[node].height = -1;
        freeList = node;
        nodeCount--;
    }

    void Tree_base::insertParticle(unsigned int particle, const std::array<double,3>& position, double radius)
    {
        // Make sure the particle doesn't already exist.
        if (particleMap.count(particle) != 0)
        {
            throw std::invalid_argument("[ERROR]: Particle already exists in tree!");
        }

        // Validate the dimensionality of the position vector.
        if (position.size() != dimension)
        {
            throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        }

        // Allocate a new node for the particle.
        unsigned int node = allocateNode();

        // AABB_base size in each dimension.
        std::array<double,3> size{};

        // Compute the AABB_base limits.
        for (unsigned int i=0;i<dimension;i++)
        {
            nodes[node].aabb_base.lowerBound[i] = position[i] - radius;
            nodes[node].aabb_base.upperBound[i] = position[i] + radius;
            size[i] = nodes[node].aabb_base.upperBound[i] - nodes[node].aabb_base.lowerBound[i];
        }

        // Fatten the AABB_base.
        for (unsigned int i=0;i<dimension;i++)
        {
            nodes[node].aabb_base.lowerBound[i] -= skinThickness * size[i];
            nodes[node].aabb_base.upperBound[i] += skinThickness * size[i];
        }
        nodes[node].aabb_base.surfaceArea = nodes[node].aabb_base.computeSurfaceArea();
        nodes[node].aabb_base.centre = nodes[node].aabb_base.computeCentre();

        // Zero the height.
        nodes[node].height = 0;

        // Insert a new leaf into the tree.
        insertLeaf(node);

        // Add the new particle to the map.
        particleMap.insert(std::unordered_map<unsigned int, unsigned int>::value_type(particle, node));

        // Store the particle index.
        nodes[node].particle = particle;
    }

    void Tree_base::insertParticle(unsigned int particle, const std::array<double,3>& lowerBound, const std::array<double,3>& upperBound)
    {
        // Make sure the particle doesn't already exist.
        if (particleMap.count(particle) != 0)
        {
            throw std::invalid_argument("[ERROR]: Particle already exists in tree!");
        }

        // Validate the dimensionality of the bounds vectors.
        if ((lowerBound.size() != dimension) || (upperBound.size() != dimension))
        {
            throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        }

        // Allocate a new node for the particle.
        unsigned int node = allocateNode();

        // AABB_base size in each dimension.
        std::array<double,3> size{};

        // Compute the AABB_base limits.
        for (unsigned int i=0;i<dimension;i++)
        {
            // Validate the bound.
            if (lowerBound[i] > upperBound[i])
            {
                throw std::invalid_argument("[ERROR]: AABB_base lower bound is greater than the upper bound!");
            }

            nodes[node].aabb_base.lowerBound[i] = lowerBound[i];
            nodes[node].aabb_base.upperBound[i] = upperBound[i];
            size[i] = upperBound[i] - lowerBound[i];
        }

        // Fatten the AABB_base.
        for (unsigned int i=0;i<dimension;i++)
        {
            nodes[node].aabb_base.lowerBound[i] -= skinThickness * size[i];
            nodes[node].aabb_base.upperBound[i] += skinThickness * size[i];
        }
        nodes[node].aabb_base.surfaceArea = nodes[node].aabb_base.computeSurfaceArea();
        nodes[node].aabb_base.centre = nodes[node].aabb_base.computeCentre();

        // Zero the height.
        nodes[node].height = 0;

        // Insert a new leaf into the tree.
        insertLeaf(node);

        // Add the new particle to the map.
        particleMap.insert(std::unordered_map<unsigned int, unsigned int>::value_type(particle, node));

        // Store the particle index.
        nodes[node].particle = particle;
    }

    unsigned int Tree_base::nParticles()
    {
        return particleMap.size();
    }

    void Tree_base::removeParticle(unsigned int particle)
    {
        // Map iterator.
        std::unordered_map<unsigned int, unsigned int>::iterator it;

        // Find the particle.
        it = particleMap.find(particle);

        // The particle doesn't exist.
        if (it == particleMap.end())
        {
            throw std::invalid_argument("[ERROR]: Invalid particle index!");
        }

        // Extract the node index.
        unsigned int node = it->second;

        // Erase the particle from the map.
        particleMap.erase(it);

        assert(node < nodeCapacity);
        assert(nodes[node].isLeaf());

        removeLeaf(node);
        freeNode(node);
    }

    void Tree_base::removeAll()
    {
        // Iterator pointing to the start of the particle map.
        std::unordered_map<unsigned int, unsigned int>::iterator it = particleMap.begin();

        // Iterate over the map.
        while (it != particleMap.end())
        {
            // Extract the node index.
            unsigned int node = it->second;

            assert(node < nodeCapacity);
            assert(nodes[node].isLeaf());

            removeLeaf(node);
            freeNode(node);

            it++;
        }

        // Clear the particle map.
        particleMap.clear();
    }

    bool Tree_base::updateParticle(unsigned int particle, std::array<double,3>& position, double radius,
                              bool alwaysReinsert)
    {
        // Validate the dimensionality of the position vector.
        if (position.size() != dimension)
        {
            throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        }

        // AABB_base bounds vectors.
        std::array<double,3> lowerBound{};
        std::array<double,3> upperBound{};

        // Compute the AABB_base limits.
        for (unsigned int i=0;i<dimension;i++)
        {
            lowerBound[i] = position[i] - radius;
            upperBound[i] = position[i] + radius;
        }

        // Update the particle.
        return updateParticle(particle, lowerBound, upperBound, alwaysReinsert);
    }

    bool Tree_base::updateParticle(unsigned int particle, std::array<double,3>& lowerBound,
                              std::array<double,3>& upperBound, bool alwaysReinsert)
    {
        // Validate the dimensionality of the bounds vectors.
        if ((lowerBound.size() != dimension) && (upperBound.size() != dimension))
        {
            throw std::invalid_argument("[ERROR]: Dimensionality mismatch!");
        }

        // Map iterator.
        std::unordered_map<unsigned int, unsigned int>::iterator it;

        // Find the particle.
        it = particleMap.find(particle);

        // The particle doesn't exist.
        if (it == particleMap.end())
        {
            throw std::invalid_argument("[ERROR]: Invalid particle index!");
        }

        // Extract the node index.
        unsigned int node = it->second;

        assert(node < nodeCapacity);
        assert(nodes[node].isLeaf());

        // AABB_base size in each dimension.
        std::array<double,3> size{};

        // Compute the AABB_base limits.
        for (unsigned int i=0;i<dimension;i++)
        {
            // Validate the bound.
            if (lowerBound[i] > upperBound[i])
            {
                throw std::invalid_argument("[ERROR]: AABB_base lower bound is greater than the upper bound!");
            }

            size[i] = upperBound[i] - lowerBound[i];
        }

        // Create the new AABB_base.
        AABB_base aabb_base(lowerBound, upperBound);

        // No need to update if the particle is still within its fattened AABB_base.
        if (!alwaysReinsert && nodes[node].aabb_base.contains(aabb_base)) return false;

        // Remove the current leaf.
        removeLeaf(node);

        // Fatten the new AABB_base.
        for (unsigned int i=0;i<dimension;i++)
        {
            aabb_base.lowerBound[i] -= skinThickness * size[i];
            aabb_base.upperBound[i] += skinThickness * size[i];
        }

        // Assign the new AABB_base.
        nodes[node].aabb_base = aabb_base;

        // Update the surface area and centroid.
        nodes[node].aabb_base.surfaceArea = nodes[node].aabb_base.computeSurfaceArea();
        nodes[node].aabb_base.centre = nodes[node].aabb_base.computeCentre();

        // Insert a new leaf node.
        insertLeaf(node);

        return true;
    }

    std::vector<unsigned int> Tree_base::query(unsigned int particle)
    {
        // Make sure that this is a valid particle.
        if (particleMap.count(particle) == 0)
        {
            throw std::invalid_argument("[ERROR]: Invalid particle index!");
        }

        // Test overlap of particle AABB_base against all other particles.
        return query(particle, nodes[particleMap.find(particle)->second].aabb_base);
    }

    std::vector<unsigned int> Tree_base::query(unsigned int particle, const AABB_base& aabb_base)
    {
        std::vector<unsigned int> stack;
        stack.reserve(256);
        stack.push_back(root);

        std::vector<unsigned int> particles;

        while (stack.size() > 0)
        {
            unsigned int node = stack.back();
            stack.pop_back();

            // Copy the AABB_base.
            AABB_base nodeAABB = nodes[node].aabb_base;

            if (node == NULL_NODE) continue;

            if (isPeriodic)
            {
                std::array<double,3> separation{};
                std::array<double,3> shift{};
                for (unsigned int i=0;i<dimension;i++)
                    separation[i] = nodeAABB.centre[i] - aabb_base.centre[i];

                bool isShifted = minimumImage(separation, shift);

                // Shift the AABB_base.
                if (isShifted)
                {
                    for (unsigned int i=0;i<dimension;i++)
                    {
                        nodeAABB.lowerBound[i] += shift[i];
                        nodeAABB.upperBound[i] += shift[i];
                    }
                }
            }

            // Test for overlap between the AABBs.
            if (aabb_base.overlaps(nodeAABB, touchIsOverlap))
            {
                // Check that we're at a leaf node.
                if (nodes[node].isLeaf())
                {
                    // Can't interact with itself.
                    if (nodes[node].particle != particle)
                    {
                        particles.push_back(nodes[node].particle);
                    }
                }
                else
                {
                    stack.push_back(nodes[node].left);
                    stack.push_back(nodes[node].right);
                }
            }
        }

        return particles;
    }

    std::vector<unsigned int> Tree_base::query(const AABB_base& aabb_base)
    {
        // Make sure the tree isn't empty.
        if (particleMap.size() == 0)
        {
            return std::vector<unsigned int>();
        }

        // Test overlap of AABB_base against all particles.
        return query(std::numeric_limits<unsigned int>::max(), aabb_base);
    }

    const AABB_base& Tree_base::getAABB(unsigned int particle)
    {
        return nodes[particleMap[particle]].aabb_base;
    }

    void Tree_base::insertLeaf(unsigned int leaf)
    {
        if (root == NULL_NODE)
        {
            root = leaf;
            nodes[root].parent = NULL_NODE;
            return;
        }

        // Find the best sibling for the node.

        AABB_base leafAABB = nodes[leaf].aabb_base;
        unsigned int index = root;

        while (!nodes[index].isLeaf())
        {
            // Extract the children of the node.
            unsigned int left  = nodes[index].left;
            unsigned int right = nodes[index].right;

            double surfaceArea = nodes[index].aabb_base.getSurfaceArea();

            AABB_base combinedAABB;
            combinedAABB.merge(nodes[index].aabb_base, leafAABB);
            double combinedSurfaceArea = combinedAABB.getSurfaceArea();

            // Cost of creating a new parent for this node and the new leaf.
            double cost = 2.0 * combinedSurfaceArea;

            // Minimum cost of pushing the leaf further down the tree.
            double inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

            // Cost of descending to the left.
            double costLeft;
            if (nodes[left].isLeaf())
            {
                AABB_base aabb_base;
                aabb_base.merge(leafAABB, nodes[left].aabb_base);
                costLeft = aabb_base.getSurfaceArea() + inheritanceCost;
            }
            else
            {
                AABB_base aabb_base;
                aabb_base.merge(leafAABB, nodes[left].aabb_base);
                double oldArea = nodes[left].aabb_base.getSurfaceArea();
                double newArea = aabb_base.getSurfaceArea();
                costLeft = (newArea - oldArea) + inheritanceCost;
            }

            // Cost of descending to the right.
            double costRight;
            if (nodes[right].isLeaf())
            {
                AABB_base aabb_base;
                aabb_base.merge(leafAABB, nodes[right].aabb_base);
                costRight = aabb_base.getSurfaceArea() + inheritanceCost;
            }
            else
            {
                AABB_base aabb_base;
                aabb_base.merge(leafAABB, nodes[right].aabb_base);
                double oldArea = nodes[right].aabb_base.getSurfaceArea();
                double newArea = aabb_base.getSurfaceArea();
                costRight = (newArea - oldArea) + inheritanceCost;
            }

            // Descend according to the minimum cost.
            if ((cost < costLeft) && (cost < costRight)) break;

            // Descend.
            if (costLeft < costRight) index = left;
            else                      index = right;
        }

        unsigned int sibling = index;

        // Create a new parent.
        unsigned int oldParent = nodes[sibling].parent;
        unsigned int newParent = allocateNode();
        nodes[newParent].parent = oldParent;
        nodes[newParent].aabb_base.merge(leafAABB, nodes[sibling].aabb_base);
        nodes[newParent].height = nodes[sibling].height + 1;

        // The sibling was not the root.
        if (oldParent != NULL_NODE)
        {
            if (nodes[oldParent].left == sibling) nodes[oldParent].left = newParent;
            else                                  nodes[oldParent].right = newParent;

            nodes[newParent].left = sibling;
            nodes[newParent].right = leaf;
            nodes[sibling].parent = newParent;
            nodes[leaf].parent = newParent;
        }
        // The sibling was the root.
        else
        {
            nodes[newParent].left = sibling;
            nodes[newParent].right = leaf;
            nodes[sibling].parent = newParent;
            nodes[leaf].parent = newParent;
            root = newParent;
        }

        // Walk back up the tree fixing heights and AABBs.
        index = nodes[leaf].parent;
        while (index != NULL_NODE)
        {
            index = balance(index);

            unsigned int left = nodes[index].left;
            unsigned int right = nodes[index].right;

            assert(left != NULL_NODE);
            assert(right != NULL_NODE);

            nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);
            nodes[index].aabb_base.merge(nodes[left].aabb_base, nodes[right].aabb_base);

            index = nodes[index].parent;
        }
    }

    void Tree_base::removeLeaf(unsigned int leaf)
    {
        if (leaf == root)
        {
            root = NULL_NODE;
            return;
        }

        unsigned int parent = nodes[leaf].parent;
        unsigned int grandParent = nodes[parent].parent;
        unsigned int sibling;

        if (nodes[parent].left == leaf) sibling = nodes[parent].right;
        else                            sibling = nodes[parent].left;

        // Destroy the parent and connect the sibling to the grandparent.
        if (grandParent != NULL_NODE)
        {
            if (nodes[grandParent].left == parent) nodes[grandParent].left = sibling;
            else                                   nodes[grandParent].right = sibling;

            nodes[sibling].parent = grandParent;
            freeNode(parent);

            // Adjust ancestor bounds.
            unsigned int index = grandParent;
            while (index != NULL_NODE)
            {
                index = balance(index);

                unsigned int left = nodes[index].left;
                unsigned int right = nodes[index].right;

                nodes[index].aabb_base.merge(nodes[left].aabb_base, nodes[right].aabb_base);
                nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);

                index = nodes[index].parent;
            }
        }
        else
        {
            root = sibling;
            nodes[sibling].parent = NULL_NODE;
            freeNode(parent);
        }
    }

    unsigned int Tree_base::balance(unsigned int node)
    {
        assert(node != NULL_NODE);

        if (nodes[node].isLeaf() || (nodes[node].height < 2))
            return node;

        unsigned int left = nodes[node].left;
        unsigned int right = nodes[node].right;

        assert(left < nodeCapacity);
        assert(right < nodeCapacity);

        int currentBalance = nodes[right].height - nodes[left].height;

        // Rotate right branch up.
        if (currentBalance > 1)
        {
            unsigned int rightLeft = nodes[right].left;
            unsigned int rightRight = nodes[right].right;

            assert(rightLeft < nodeCapacity);
            assert(rightRight < nodeCapacity);

            // Swap node and its right-hand child.
            nodes[right].left = node;
            nodes[right].parent = nodes[node].parent;
            nodes[node].parent = right;

            // The node's old parent should now point to its right-hand child.
            if (nodes[right].parent != NULL_NODE)
            {
                if (nodes[nodes[right].parent].left == node) nodes[nodes[right].parent].left = right;
                else
                {
                    assert(nodes[nodes[right].parent].right == node);
                    nodes[nodes[right].parent].right = right;
                }
            }
            else root = right;

            // Rotate.
            if (nodes[rightLeft].height > nodes[rightRight].height)
            {
                nodes[right].right = rightLeft;
                nodes[node].right = rightRight;
                nodes[rightRight].parent = node;
                nodes[node].aabb_base.merge(nodes[left].aabb_base, nodes[rightRight].aabb_base);
                nodes[right].aabb_base.merge(nodes[node].aabb_base, nodes[rightLeft].aabb_base);

                nodes[node].height = 1 + std::max(nodes[left].height, nodes[rightRight].height);
                nodes[right].height = 1 + std::max(nodes[node].height, nodes[rightLeft].height);
            }
            else
            {
                nodes[right].right = rightRight;
                nodes[node].right = rightLeft;
                nodes[rightLeft].parent = node;
                nodes[node].aabb_base.merge(nodes[left].aabb_base, nodes[rightLeft].aabb_base);
                nodes[right].aabb_base.merge(nodes[node].aabb_base, nodes[rightRight].aabb_base);

                nodes[node].height = 1 + std::max(nodes[left].height, nodes[rightLeft].height);
                nodes[right].height = 1 + std::max(nodes[node].height, nodes[rightRight].height);
            }

            return right;
        }

        // Rotate left branch up.
        if (currentBalance < -1)
        {
            unsigned int leftLeft = nodes[left].left;
            unsigned int leftRight = nodes[left].right;

            assert(leftLeft < nodeCapacity);
            assert(leftRight < nodeCapacity);

            // Swap node and its left-hand child.
            nodes[left].left = node;
            nodes[left].parent = nodes[node].parent;
            nodes[node].parent = left;

            // The node's old parent should now point to its left-hand child.
            if (nodes[left].parent != NULL_NODE)
            {
                if (nodes[nodes[left].parent].left == node) nodes[nodes[left].parent].left = left;
                else
                {
                    assert(nodes[nodes[left].parent].right == node);
                    nodes[nodes[left].parent].right = left;
                }
            }
            else root = left;

            // Rotate.
            if (nodes[leftLeft].height > nodes[leftRight].height)
            {
                nodes[left].right = leftLeft;
                nodes[node].left = leftRight;
                nodes[leftRight].parent = node;
                nodes[node].aabb_base.merge(nodes[right].aabb_base, nodes[leftRight].aabb_base);
                nodes[left].aabb_base.merge(nodes[node].aabb_base, nodes[leftLeft].aabb_base);

                nodes[node].height = 1 + std::max(nodes[right].height, nodes[leftRight].height);
                nodes[left].height = 1 + std::max(nodes[node].height, nodes[leftLeft].height);
            }
            else
            {
                nodes[left].right = leftRight;
                nodes[node].left = leftLeft;
                nodes[leftLeft].parent = node;
                nodes[node].aabb_base.merge(nodes[right].aabb_base, nodes[leftLeft].aabb_base);
                nodes[left].aabb_base.merge(nodes[node].aabb_base, nodes[leftRight].aabb_base);

                nodes[node].height = 1 + std::max(nodes[right].height, nodes[leftLeft].height);
                nodes[left].height = 1 + std::max(nodes[node].height, nodes[leftRight].height);
            }

            return left;
        }

        return node;
    }

    unsigned int Tree_base::computeHeight() const
    {
        return computeHeight(root);
    }

    unsigned int Tree_base::computeHeight(unsigned int node) const
    {
        assert(node < nodeCapacity);

        if (nodes[node].isLeaf()) return 0;

        unsigned int height1 = computeHeight(nodes[node].left);
        unsigned int height2 = computeHeight(nodes[node].right);

        return 1 + std::max(height1, height2);
    }

    unsigned int Tree_base::getHeight() const
    {
        if (root == NULL_NODE) return 0;
        return nodes[root].height;
    }

    unsigned int Tree_base::getNodeCount() const
    {
        return nodeCount;
    }

    unsigned int Tree_base::computeMaximumBalance() const
    {
        unsigned int maxBalance = 0;
        for (unsigned int i=0; i<nodeCapacity; i++)
        {
            if (nodes[i].height <= 1)
                continue;

            assert(nodes[i].isLeaf() == false);

            unsigned int balance = std::abs(nodes[nodes[i].left].height - nodes[nodes[i].right].height);
            maxBalance = std::max(maxBalance, balance);
        }

        return maxBalance;
    }

    double Tree_base::computeSurfaceAreaRatio() const
    {
        if (root == NULL_NODE) return 0.0;

        double rootArea = nodes[root].aabb_base.computeSurfaceArea();
        double totalArea = 0.0;

        for (unsigned int i=0; i<nodeCapacity;i++)
        {
            if (nodes[i].height < 0) continue;

            totalArea += nodes[i].aabb_base.computeSurfaceArea();
        }

        return totalArea / rootArea;
    }

    void Tree_base::validate() const
    {
#ifndef NDEBUG
        validateStructure(root);
        validateMetrics(root);

        unsigned int freeCount = 0;
        unsigned int freeIndex = freeList;

        while (freeIndex != NULL_NODE)
        {
            assert(freeIndex < nodeCapacity);
            freeIndex = nodes[freeIndex].next;
            freeCount++;
        }

        assert(getHeight() == computeHeight());
        assert((nodeCount + freeCount) == nodeCapacity);
#endif
    }

    void Tree_base::rebuild()
    {
        std::vector<unsigned int> nodeIndices(nodeCount);
        unsigned int count = 0;

        for (unsigned int i=0;i<nodeCapacity;i++)
        {
            // Free node.
            if (nodes[i].height < 0) continue;

            if (nodes[i].isLeaf())
            {
                nodes[i].parent = NULL_NODE;
                nodeIndices[count] = i;
                count++;
            }
            else freeNode(i);
        }

        while (count > 1)
        {
            double minCost = std::numeric_limits<double>::max();
            int iMin = -1, jMin = -1;

            for (unsigned int i=0;i<count;i++)
            {
                AABB_base aabbi = nodes[nodeIndices[i]].aabb_base;

                for (unsigned int j=i+1;j<count;j++)
                {
                    AABB_base aabbj = nodes[nodeIndices[j]].aabb_base;
                    AABB_base aabb_base;
                    aabb_base.merge(aabbi, aabbj);
                    double cost = aabb_base.getSurfaceArea();

                    if (cost < minCost)
                    {
                        iMin = i;
                        jMin = j;
                        minCost = cost;
                    }
                }
            }

            unsigned int index1 = nodeIndices[iMin];
            unsigned int index2 = nodeIndices[jMin];

            unsigned int parent = allocateNode();
            nodes[parent].left = index1;
            nodes[parent].right = index2;
            nodes[parent].height = 1 + std::max(nodes[index1].height, nodes[index2].height);
            nodes[parent].aabb_base.merge(nodes[index1].aabb_base, nodes[index2].aabb_base);
            nodes[parent].parent = NULL_NODE;

            nodes[index1].parent = parent;
            nodes[index2].parent = parent;

            nodeIndices[jMin] = nodeIndices[count-1];
            nodeIndices[iMin] = parent;
            count--;
        }

        root = nodeIndices[0];

        validate();
    }

    void Tree_base::validateStructure(unsigned int node) const
    {
        if (node == NULL_NODE) return;

        if (node == root) assert(nodes[node].parent == NULL_NODE);

        unsigned int left = nodes[node].left;
        unsigned int right = nodes[node].right;

        if (nodes[node].isLeaf())
        {
            assert(left == NULL_NODE);
            assert(right == NULL_NODE);
            assert(nodes[node].height == 0);
            return;
        }

        assert(left < nodeCapacity);
        assert(right < nodeCapacity);

        assert(nodes[left].parent == node);
        assert(nodes[right].parent == node);

        validateStructure(left);
        validateStructure(right);
    }

    void Tree_base::validateMetrics(unsigned int node) const
    {
        if (node == NULL_NODE) return;

        unsigned int left = nodes[node].left;
        unsigned int right = nodes[node].right;

        if (nodes[node].isLeaf())
        {
            assert(left == NULL_NODE);
            assert(right == NULL_NODE);
            assert(nodes[node].height == 0);
            return;
        }

        assert(left < nodeCapacity);
        assert(right < nodeCapacity);

        int height1 = nodes[left].height;
        int height2 = nodes[right].height;
        int height = 1 + std::max(height1, height2);
        (void)height; // Unused variable in Release build
        assert(nodes[node].height == height);

        AABB_base aabb_base;
        aabb_base.merge(nodes[left].aabb_base, nodes[right].aabb_base);

        for (unsigned int i=0;i<dimension;i++)
        {
            assert(aabb_base.lowerBound[i] == nodes[node].aabb_base.lowerBound[i]);
            assert(aabb_base.upperBound[i] == nodes[node].aabb_base.upperBound[i]);
        }

        validateMetrics(left);
        validateMetrics(right);
    }

    void Tree_base::periodicBoundaries(std::array<double,3>& position)
    {
        for (unsigned int i=0;i<dimension;i++)
        {
            if (position[i] < 0)
            {
                position[i] += boxSize[i];
            }
            else
            {
                if (position[i] >= boxSize[i])
                {
                    position[i] -= boxSize[i];
                }
            }
        }
    }

    bool Tree_base::minimumImage(std::array<double,3>& separation, std::array<double,3>& shift)
    {
        bool isShifted = false;

        for (unsigned int i=0;i<dimension;i++)
        {
            if (separation[i] < negMinImage[i])
            {
                separation[i] += periodicity[i]*boxSize[i];
                shift[i] = periodicity[i]*boxSize[i];
                isShifted = true;
            }
            else
            {
                if (separation[i] >= posMinImage[i])
                {
                    separation[i] -= periodicity[i]*boxSize[i];
                    shift[i] = -static_cast<int>(periodicity[i])*boxSize[i];
                    isShifted = true;
                }
            }
        }

        return isShifted;
    }
}
