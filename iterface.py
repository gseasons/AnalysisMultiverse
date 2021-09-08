#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 14:10:12 2021

@author: grahamseasons
"""
from nipype.interfaces.base import (
    DynamicTraitedSpec,
    isdefined,
)
from nipype.interfaces.io import IOBase, add_traits

#change so that if just a single item -> returns it

class IdentityIterface(IOBase):
    """Basic interface class generates identity mappings

    Examples
    --------

    >>> from nipype.interfaces.utility import IdentityInterface
    >>> ii = IdentityInterface(fields=['a', 'b'], mandatory_inputs=False)
    >>> ii.inputs.a
    <undefined>

    >>> ii.inputs.a = 'foo'
    >>> out = ii._outputs()
    >>> out.a
    <undefined>

    >>> out = ii.run()
    >>> out.outputs.a
    'foo'

    >>> ii2 = IdentityInterface(fields=['a', 'b'], mandatory_inputs=True)
    >>> ii2.inputs.a = 'foo'
    >>> out = ii2.run() # doctest: +SKIP
    ValueError: IdentityInterface requires a value for input 'b' because it was listed in 'fields' Interface IdentityInterface failed to run.
    """

    input_spec = DynamicTraitedSpec
    output_spec = DynamicTraitedSpec

    def __init__(self, fields=None, index=None, mandatory_inputs=True, **inputs):
        super().__init__(**inputs)
        if fields is None or not fields:
            raise ValueError("Identity Interface fields must be a non-empty list")
        # Each input must be in the fields.
        for in_field in inputs:
            if in_field not in fields:
                raise ValueError(
                    "Identity Interface input is not in the fields: %s" % in_field
                )
        self._fields = fields
        self._index = index
        self._mandatory_inputs = mandatory_inputs
        add_traits(self.inputs, fields)#, index)
        # Adding any traits wipes out all input values set in superclass initialization,
        # even it the trait is not in the add_traits argument. The work-around is to reset
        # the values after adding the traits.
        self.inputs.trait_set(**inputs)

    def _add_output_traits(self, base):
        return add_traits(base, self._fields)
    def _list_outputs(self):
        # manual mandatory inputs check
        if self._fields and self._mandatory_inputs:
            flag = True
            if 'index' in self._fields:
                for i, key in enumerate(self._fields):
                    value = getattr(self.inputs, key)
                    if key == 'index':
                        self._index = value
                        continue
                    
                    if flag and type(value) == list:
                        old_val = len(value)
                        flag = False
                    if not isdefined(value):
                        msg = (
                            "%s requires a value for input '%s' because it was listed in 'fields'. \
                        You can turn off mandatory inputs checking by passing mandatory_inputs = False to the constructor."
                            % (self.__class__.__name__, key)
                        )
                        raise ValueError(msg)
                    if type(value) != list:
                        continue
                    if not (type(value) == list):
                        msg = (
                            "%s needs to be a list as IdentityIterface is designed to iterate through lists."
                            % (self.__class__.__name__)
                        )
                        raise ValueError(msg)
                    if old_val != len(value):
                        msg = (
                            "All inputs need to be of equal length."
                        )
                        raise ValueError(msg)
                    old_val = len(value)
        if self._index == None:
            msg = (
                "index must be specified."
                )
            raise ValueError(msg)
        if self._index >= old_val:
            msg = (
                "index out of range."
                )
            raise ValueError(msg)
            

        outputs = self._outputs().get()
        for key in self._fields:
            val = getattr(self.inputs, key)
            if isdefined(val):
                if type(val) != list:
                    outputs[key] = val
                else:
                    outputs[key] = val[self._index]
        return outputs
