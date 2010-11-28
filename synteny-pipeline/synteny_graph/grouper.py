# Helper functions in the blast filtering
def gene_name(st):
    # this is ugly, but different groups are inconsistent
    # with how the alternative splicings are named;
    # mostly it can be done by removing the suffix
    # except for papaya (evm...) and maize (somewhat complicated)
    # TODO: maize rules
    if st.startswith("ev"):
        return st
    if st.startswith("Os"):
        return st.rsplit("-",1)[0]
    return st.rsplit(".", 1)[0]

# colored output
def colored_text(st, color="yellow"):
    color_supported = "gray|red|green|yellow|blue|magenta|cyan|white|crimson"
    color_supported = color_supported.split("|")
    if color not in color_supported:
        return st

    color_instances = dict(zip(color_supported, range(30, 39)))
    # \033[1;??m, ?? can range from 30 to 48
    return '\033[1;%dm%s\033[1;m' % (color_instances[color], st)

"""
Disjoint set data structure [http://code.activestate.com/recipes/387776/]
"""
class Grouper(object):
   """This class provides a lightweight way to group arbitrary objects
together into disjoint sets when a full-blown graph data structure
would be overkill.

Objects can be joined using .join(), tested for connectedness
using .joined(), and all disjoint sets can be retreived using
.get().

The objects being joined must be hashable.

>>> g = Grouper()
>>> g.join('a', 'b')
>>> g.join('b', 'c')
>>> g.join('d', 'e')
>>> list(g)
[['a', 'b', 'c'], ['d', 'e']]
>>> g.joined('a', 'b')
True
>>> g.joined('a', 'c')
True
>>> 'f' in g
False
>>> g.joined('a', 'd')
False"""   
   def __init__(self, init=[]):
      mapping = self._mapping = {}
      for x in init:
         mapping[x] = [x]
        
   def join(self, a, *args):
      """Join given arguments into the same set. Accepts one or more arguments."""
      mapping = self._mapping
      set_a = mapping.setdefault(a, [a])

      for arg in args:
         set_b = mapping.get(arg)
         if set_b is None:
            set_a.append(arg)
            mapping[arg] = set_a
         elif set_b is not set_a:
            if len(set_b) > len(set_a):
               set_a, set_b = set_b, set_a
            set_a.extend(set_b)
            for elem in set_b:
               mapping[elem] = set_a

   def joined(self, a, b):
      """Returns True if a and b are members of the same set."""
      mapping = self._mapping
      try:
          return mapping[a] is mapping[b]
      except KeyError:
          return False

   def __iter__(self):
      """Returns an iterator returning each of the disjoint sets as a list."""
      seen = set()
      for elem, group in self._mapping.iteritems():
          if elem not in seen:
              yield group
              seen.update(group)

   def __getitem__(self, key):
       """Returns the set that a certain key belongs."""
       return tuple(self._mapping[key])

   def __contains__(self, key):
       return key in self._mapping

   def __len__(self):
       #return len(list(self))
       group = set()
       for v in self._mapping.values():
           group.update([tuple(v)])
       return len(group)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
