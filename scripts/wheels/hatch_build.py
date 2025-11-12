import os
import shutil

from pathlib import Path
from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from hatchling.metadata.plugin.interface import MetadataHookInterface

class CustomHook(BuildHookInterface):
    def initialize(self, version, build_data):
        # Generate Wheel tags. This is necessary because Kratos is platform specific.
        build_data["infer_tag"] = False
        build_data["pure_pyton"] = False

        abi_tag = os.environ['ABI_TAG']

        build_data["tag"] = f"{abi_tag}"

class CustomMeta(MetadataHookInterface):
    def update(self, metadata):
        # Fill version globaly for all projects
        metadata["version"] = os.environ["KRATOS_VERSION"]