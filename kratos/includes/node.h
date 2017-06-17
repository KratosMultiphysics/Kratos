
        
#ifdef KRATOS_DEBUG
        if(rDofVariable.Key() == 0) 
        {
            KRATOS_ERROR << "Variable  " << rDofVariable << " has key zero key when adding Dof for node " << this->Id() << std::endl;
        }
#endif
        
#ifdef KRATOS_DEBUG
        if(rDofVariable.Key() == 0) 
        {
            KRATOS_ERROR << "Variable  " << rDofVariable << " has key zero key when adding Dof for node " << this->Id() << std::endl;
        }
        if(rDofReaction.Key() == 0) 
        {
            KRATOS_ERROR << "Reaction  " << rDofReaction << " has key zero when adding reactions for node " << this->Id() << std::endl;
        }
#endif
        
#ifdef KRATOS_DEBUG
        if(rDofVariable.Key() == 0)
        {
            KRATOS_ERROR << "Variable  " << rDofVariable << " has key zero key when adding Dof for node " << this->Id() << std::endl;
        }
#endif
        typename DofsContainerType::iterator itDoF = mDofs.find(rDofVariable);
        if(itDoF != mDofs.end())
        {
            return *itDoF;
        }
            
        typename DofType::Pointer pNewDoF =  boost::make_shared<DofType>(Id(), &mSolutionStepsNodalData, rDofVariable);
        mDofs.insert(mDofs.begin(), pNewDoF);
        
#ifdef KRATOS_DEBUG
        if(rDofVariable.Key() == 0) 
        {
            KRATOS_ERROR << "Variable  " << rDofVariable << " has key zero key when adding Dof for node " << this->Id() << std::endl;
        }
        if(rDofReaction.Key() == 0) 
        {
            KRATOS_ERROR << "Reaction  " << rDofReaction << " has key zero when adding reactions for node " << this->Id() << std::endl;
        }
#endif
        typename DofsContainerType::iterator itDoF = mDofs.find(rDofVariable);
        if(itDoF != mDofs.end())
        {
            itDoF->SetReaction(rDofReaction);
            return *itDoF;
        }