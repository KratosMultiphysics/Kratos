---
title: Check Skin
keywords: process core
tags: [check skin process]
sidebar: kratos_core_processes
summary: 
---

# Check Skin

## Description

**Warning**: This process does not follow the standard Process interface.

This function verifies that the skin has no holes nor overlapped geometries this is accomplished by storing all of the edges in the model in a hash map and verifying that no edge appears more than twice (which would imply an overlap) or less than once (which would imply a gap in the skin)

Please note that checks process do not produce an output, and will only print its result.

## Parameters & Defaults

This process is not configurable