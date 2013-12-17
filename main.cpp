#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <vector>

// ---- TomF's vertex cache optimizer
// my old Implementation from Werkkzeug4/Altona copy & pasted

struct VCacheVert
{
  int CachePos;      // its position in the cache (-1 if not in)
  int Score;         // its score (higher=better)
  int TrisLeft;      // # of not-yet-used tris
  int *TriList;      // list of triangle indices
  int OpenPos;       // position in "open vertex" list
};

struct VCacheTri
{
  int Score;         // current score (-1 if already done)
  int Inds[3];       // vertex indices
};

void OptimizeIndexOrder(int *IndexBuffer,int IndexCount,int VertexCount)
{
  VCacheVert *verts = new VCacheVert[VertexCount];
  for(int i=0;i<VertexCount;i++)
  {
    verts[i].CachePos = -1;
    verts[i].Score = 0;
    verts[i].TrisLeft = 0;
    verts[i].TriList = 0;
    verts[i].OpenPos = -1;
  }

  // prepare triangles
  int nTris = IndexCount/3;
  VCacheTri *tris = new VCacheTri[nTris];
  int *indPtr = IndexBuffer;

  for(int i=0;i<nTris;i++)
  {
    tris[i].Score = 0;

    for(int j=0;j<3;j++)
    {
      int ind = *indPtr++;
      tris[i].Inds[j] = ind;
      verts[ind].TrisLeft++;
    }
  }

  // alloc space for vert->tri indices
  int *vertTriInd = new int[nTris*3];
  int *vertTriPtr = vertTriInd;

  for(int i=0;i<VertexCount;i++)
  {
    verts[i].TriList = vertTriPtr;
    vertTriPtr += verts[i].TrisLeft;
    verts[i].TrisLeft = 0;
  }

  // make vert->tri tables
  for(int i=0;i<nTris;i++)
  {
    for(int j=0;j<3;j++)
    {
      int ind = tris[i].Inds[j];
      verts[ind].TriList[verts[ind].TrisLeft] = i;
      verts[ind].TrisLeft++;
    }
  }

  // open vertices
  int *openVerts = new int[VertexCount];
  int openCount = 0;

  // the cache
  static const int cacheSize = 32;
  static const int maxValence = 15;
  int cache[cacheSize+3];
  int pos2Score[cacheSize];
  int val2Score[maxValence+1];

  for(int i=0;i<cacheSize+3;i++)
    cache[i] = -1;

  for(int i=0;i<cacheSize;i++)
  {
    float score = (i<3) ? 0.75f : powf(1.0f - (i-3)/float(cacheSize-3),1.5f);
    pos2Score[i] = (int) (score * 65536.0f + 0.5f);
  }

  val2Score[0] = 0;
  for(int i=1;i<16;i++)
  {
    float score = 2.0f / sqrtf((float) i);
    val2Score[i] = (int) (score * 65536.0f + 0.5f);
  }

  // outer loop: find triangle to start with
  indPtr = IndexBuffer;
  int seedPos = 0;

  while(1)
  {
    int seedScore = -1;
    int seedTri = -1;

    // if there are open vertices, search them for the seed triangle
    // which maximum score.
    for(int i=0;i<openCount;i++)
    {
      VCacheVert *vert = &verts[openVerts[i]];

      for(int j=0;j<vert->TrisLeft;j++)
      {
        int triInd = vert->TriList[j];
        VCacheTri *tri = &tris[triInd];

        if(tri->Score > seedScore)
        {
          seedScore = tri->Score;
          seedTri = triInd;
        }
      }
    }

    // if we haven't found a seed triangle yet, there are no open
    // vertices and we can pick any triangle
    if(seedTri == -1)
    {
      while(seedPos < nTris && tris[seedPos].Score<0)
        seedPos++;

      if(seedPos == nTris) // no triangle left, we're done!
        break;

      seedTri = seedPos;
    }

    // the main loop.
    int bestTriInd = seedTri;
    while(bestTriInd != -1)
    {
      VCacheTri *bestTri = &tris[bestTriInd];

      // mark this triangle as used, remove it from the "remaining tris"
      // list of the vertices it uses, and add it to the index buffer.
      bestTri->Score = -1;

      for(int j=0;j<3;j++)
      {
        int vertInd = bestTri->Inds[j];
        *indPtr++ = vertInd;

        VCacheVert *vert = &verts[vertInd];
        
        // find this triangles' entry
        int k = 0;
        while(vert->TriList[k] != bestTriInd)
        {
          assert(k < vert->TrisLeft);
          k++;
        }

        // swap it to the end and decrement # of tris left
        if(--vert->TrisLeft)
          std::swap(vert->TriList[k],vert->TriList[vert->TrisLeft]);
        else if(vert->OpenPos >= 0)
          std::swap(openVerts[vert->OpenPos],openVerts[--openCount]);
      }

      // update cache status
      cache[cacheSize] = cache[cacheSize+1] = cache[cacheSize+2] = -1;

      for(int j=0;j<3;j++)
      {
        int ind = bestTri->Inds[j];
        cache[cacheSize+2] = ind;

        // find vertex index
        int pos;
        for(pos=0;cache[pos]!=ind;pos++);

        // move to front
        for(int k=pos;k>0;k--)
          cache[k] = cache[k-1];

        cache[0] = ind;

        // remove sentinel if it wasn't used
        if(pos!=cacheSize+2)
          cache[cacheSize+2] = -1;
      }

      // update vertex scores
      for(int i=0;i<cacheSize+3;i++)
      {
        int vertInd = cache[i];
        if(vertInd == -1)
          continue;

        VCacheVert *vert = &verts[vertInd];

        vert->Score = val2Score[std::min(vert->TrisLeft,maxValence)];
        if(i < cacheSize)
        {
          vert->CachePos = i;
          vert->Score += pos2Score[i];
        }
        else
          vert->CachePos = -1;

        // also add to open vertices list if the vertex is indeed open
        if(vert->OpenPos<0 && vert->TrisLeft)
        {
          vert->OpenPos = openCount;
          openVerts[openCount++] = vertInd;
        }
      }

      // update triangle scores, find new best triangle
      int bestTriScore = -1;
      bestTriInd = -1;

      for(int i=0;i<cacheSize;i++)
      {
        if(cache[i] == -1)
          continue;

        const VCacheVert *vert = &verts[cache[i]];

        for(int j=0;j<vert->TrisLeft;j++)
        {
          int triInd = vert->TriList[j];
          VCacheTri *tri = &tris[triInd];

          assert(tri->Score != -1);

          int score = 0;
          for(int k=0;k<3;k++)
            score += verts[tri->Inds[k]].Score;

          tri->Score = score;
          if(score > bestTriScore)
          {
            bestTriScore = score;
            bestTriInd = triInd;
          }
        }
      }
    }
  }

  // cleanup
  delete[] verts;
  delete[] tris;
  delete[] vertTriInd;
  delete[] openVerts;
}

/****************************************************************************/

void DumpCacheEfficiency(const int *inds,int indCount)
{
  static const int cacheSize = 16;
  static const bool isFIFO = true;

  int misses = 0;
  int cache[cacheSize+1];
  int wrPos = 0;

  if(!indCount)
    return;

  // initialize cache (we simulate a FIFO here)
  for(int i=0;i<cacheSize;i++)
    cache[i] = -1;

  // simulate
  for(int i=0;i<indCount;i++)
  {
    int ind = inds[i];
    cache[cacheSize] = ind;

    // find in cache
    int cachePos;
    for(cachePos=0;cache[cachePos] != inds[i];cachePos++);
    misses += cachePos == cacheSize;

    if(isFIFO)
    {
      cache[wrPos] = ind;
      if(++wrPos == cacheSize)
        wrPos = 0;
    }
    else
    {
      // move to front
      for(int j=cachePos;j>0;j--)
        cache[j] = cache[j-1];

      cache[0] = ind;
    }
  }
    
  // print results
  float ACMR = misses * 3.0f / indCount;
  printf("ACMR: %.3f (%d-entry %s)\n", ACMR, cacheSize, isFIFO ? "FIFO" : "LRU");
}

// ---- The actual mesh processing code

struct Vert
{
    float x, y, z;
};

struct Mesh
{
    std::vector<Vert> verts;
    std::vector<int> inds;
};

static void renumber_verts(Mesh *m)
{
    std::vector<int> old_to_new;
    std::vector<Vert> new_verts;

    old_to_new.resize(m->verts.size(), -1);
    new_verts.reserve(m->verts.size());

    // construct mapping and new VB
    for (int &ind : m->inds)
    {
        if (old_to_new[ind] == -1) // first time we see this vert
        {
            old_to_new[ind] = (int)new_verts.size();
            new_verts.push_back(m->verts[ind]);
        }
        ind = old_to_new[ind];
    }

    // switch to new VB
    m->verts.swap(new_verts);
}

static void read_mesh(Mesh *mesh, const char *filename)
{
    mesh->verts.clear();
    mesh->inds.clear();

    FILE *f = fopen(filename, "rb");
    int nverts, ninds;

    fread(&nverts, sizeof(int), 1, f);
    fread(&ninds, sizeof(int), 1, f);

    mesh->verts.resize(nverts);
    mesh->inds.resize(ninds);

    fread(&mesh->verts[0], sizeof(Vert), nverts, f);
    fread(&mesh->inds[0], sizeof(int), ninds, f);

    fclose(f);
}

static void write_mesh(const char *filename, const std::vector<Vert>& verts, const std::vector<int>& inds)
{
    FILE *f = fopen(filename, "wb");

    int nverts = verts.size();
    int ninds = inds.size();
    fwrite(&nverts, sizeof(int), 1, f);
    fwrite(&ninds, sizeof(int), 1, f);

    fwrite(&verts[0], sizeof(Vert), nverts, f);
    fwrite(&inds[0], sizeof(int), ninds, f);

    fclose(f);
}

static bool tri_has_edge(const int *tri, int a, int b)
{
    return tri[0] == a && tri[1] == b ||
        tri[1] == a && tri[2] == b ||
        tri[2] == a && tri[0] == b;
}

static void analyze_inds(const std::vector<int>& inds)
{
    int num_single = 0;
    int num_pairs = 0;

    for (size_t base = 0; base < inds.size(); )
    {
        bool is_single = true;
        if (base + 3 < inds.size()) // we can look at next tri
        {
            const int *tri = &inds[base];
            const int *next = &inds[base + 3];

            // does the next tri contain the opposite
            // of at least one of our edges?
            if (tri_has_edge(next, tri[1], tri[0]) ||
                tri_has_edge(next, tri[2], tri[1]) ||
                tri_has_edge(next, tri[0], tri[2]))
                is_single = false;
        }

        if (is_single)
        {
            num_single++;
            base += 3;
        }
        else
        {
            num_pairs++;
            base += 6;
        }
    }

    int count = inds.size();
    int fancy_ib_size = 3*num_single + 4*num_pairs;

    printf("%d paired tris, %d single\n", num_pairs * 2, num_single);
    printf("IB inds: list=%d, fancy=%d (%.2f%%)\n", count, fancy_ib_size, 100.0f * (fancy_ib_size - count) / count);
}

static void add_tri(std::vector<int>& out_inds, int a, int b, int c)
{
    assert(a >= b);
    out_inds.push_back(a);
    out_inds.push_back(b);
    out_inds.push_back(c);
}

//  c__
//  |  ^^--b
//   |    / \
//    |  /   \
//     |/   __d
//     a--^^
//
// (a,b,c,d) (which unpacks to a,b,c, a,d,b)
static void add_double_tri(std::vector<int>& out_inds, int a, int b, int c, int d)
{
    assert(a < b);
    out_inds.push_back(a);
    out_inds.push_back(b);
    out_inds.push_back(c);
    out_inds.push_back(d);
}

static bool try_merge_with_next(std::vector<int>& out_inds, const std::vector<int>& inds, size_t base)
{
    // is there even a next tri?
    if (base + 3 >= inds.size())
        return false;

    // is this tri degenerate?
    const int *tri = &inds[base];
    if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0])
        return false;

    // does the next tri contain the opposite of at least one
    // of our edges?
    const int *next = &inds[base + 3];

    // go through 3 edges of tri
    for (int i = 0; i < 3; i++)
    {
        // try to find opposite of edge ab, namely ba.
        int a = tri[i];
        int b = tri[(i + 1) % 3];
        int c = tri[(i + 2) % 3];

        for (int j = 0; j < 3; j++)
        {
            if (next[j] == b && next[(j + 1) % 3] == a)
            {
                int d = next[(j + 2) % 3];

                if (a < b)
                    add_double_tri(out_inds, a, b, c, d);
                else // must be c > a, since we checked that a != c above; this ends up swapping two tris.
                    add_double_tri(out_inds, b, a, d, c);

                return true;
            }
        }
    }

    return false;
}

static void pack_inds(std::vector<int>& out_inds, const std::vector<int>& inds)
{
    for (size_t base = 0; base < inds.size(); )
    {
        if (try_merge_with_next(out_inds, inds, base))
            base += 6; // consume 2 tris
        else
        {
            const int *tri = &inds[base];
            if (tri[0] >= tri[1])
                add_tri(out_inds, tri[0], tri[1], tri[2]);
            else if (tri[1] >= tri[2])
                add_tri(out_inds, tri[1], tri[2], tri[0]);
            else
            {
                // must have tri[2] >= tri[0],
                // otherwise we'd have tri[0] < tri[1] < tri[2] < tri[0] (contradiction)
                add_tri(out_inds, tri[2], tri[0], tri[1]);
            }

            base += 3;
        }
    }
}

static void unpack_inds(std::vector<int>& out_inds, const std::vector<int>& inds)
{
    for (size_t base = 0; base < inds.size(); )
    {
        int a = inds[base++];
        int b = inds[base++];
        int c = inds[base++];
        out_inds.push_back(a); out_inds.push_back(b); out_inds.push_back(c);

        if (a < b) // two tris: (a,b,c), (a,d,b)
        {
            int d = inds[base++];
            out_inds.push_back(a); out_inds.push_back(d); out_inds.push_back(b);
        }
    }
}

static bool tris_same(const int* tri0, const int* tri1)
{
    // stupid way to check this, but it's simple
    int a0 = tri0[0], b0 = tri0[1], c0 = tri0[2];
    int a1 = tri1[0], b1 = tri1[1], c1 = tri1[2];

    return
        (a0 == a1 && b0 == b1 && c0 == c1) ||
        (a0 == b1 && b0 == c1 && c0 == a1) ||
        (a0 == c1 && b0 == a1 && c0 == b1);
}

// this checks that the index buffers are the same, except for the permutations
// we expect from the reordering.
static bool inds_match(const std::vector<int>& inds0, const std::vector<int>& inds1)
{
    if (inds0.size() != inds1.size())
        return false;

    size_t count = inds0.size();

    // invariant: all tris up to "base" matched successfully
    for (size_t base = 0; base < count; )
    {
        if (tris_same(&inds0[base], &inds1[base]))
        {
            // tris at this position match. advance.
            base += 3;
        }
        else if (base + 3 < count && tris_same(&inds0[base + 3], &inds1[base]) && tris_same(&inds0[base], &inds1[base + 3]))
        {
            // pair of tris that matches but has swapped order
            base += 6;
        }
        else
            return false; // no match, something went wrong!
    }

    return true;
}

static void watermark_transform(std::vector<int>& out_inds, const std::vector<int>& in_inds, int max_step)
{
    int hi = max_step - 1; // high watermark
    out_inds.clear();
    out_inds.reserve(in_inds.size());
    for (int v : in_inds)
    {
        assert(v <= hi);
        out_inds.push_back(hi - v);
        hi = std::max(hi, v + max_step);
    }
}

static void reverse_watermark_transform(std::vector<int>& out_inds, const std::vector<int>& in_inds, int max_step)
{
    int hi = max_step - 1; // high watermark
    out_inds.clear();
    out_inds.reserve(in_inds.size());
    for (int v : in_inds)
    {
        assert(v <= hi);
        v = hi - v;
        out_inds.push_back(v);
        hi = std::max(hi, v + max_step);
    }
}

static void watermark_transform_and_check(std::vector<int>& out_inds, const std::vector<int>& in_inds, int max_step)
{
    std::vector<int> temp;
    watermark_transform(out_inds, in_inds, max_step);
    reverse_watermark_transform(temp, out_inds, max_step);
    assert(std::equal(temp.begin(), temp.end(), in_inds.begin()));
}

int main()
{
    Mesh m;
    read_mesh(&m, "Armadillo.bin");

    printf("%d verts, %d inds.\n", (int) m.verts.size(), (int) m.inds.size());

    printf("before:\n");
    DumpCacheEfficiency(&m.inds[0], m.inds.size());
    analyze_inds(m.inds);

    OptimizeIndexOrder(&m.inds[0], m.inds.size(), m.verts.size());
    renumber_verts(&m);
    printf("after:\n");
    DumpCacheEfficiency(&m.inds[0], m.inds.size());
    analyze_inds(m.inds);

    std::vector<int> packed_inds, unpacked_inds;
    pack_inds(packed_inds, m.inds);
    printf("%d inds packed\n", (int) packed_inds.size());

    unpack_inds(unpacked_inds, packed_inds);
    printf("%d inds unpacked\n", (int) unpacked_inds.size());

    // check that unpacking works
    if (!inds_match(m.inds, unpacked_inds))
        printf("ERROR: something went wrong!\n");
    else
        printf("index buffers match.\n");

    DumpCacheEfficiency(&unpacked_inds[0], unpacked_inds.size());

    write_mesh("Armadillo_vcache.bin", m.verts, m.inds);

    std::vector<int> watermark_inds;
    watermark_transform_and_check(watermark_inds, m.inds, 1);
    write_mesh("Armadillo_vcache_wmark.bin", m.verts, watermark_inds);

    write_mesh("Armadillo_post.bin", m.verts, packed_inds);
    
    std::vector<int> watermark_packed_inds;
    watermark_transform_and_check(watermark_packed_inds, packed_inds, 3);
    write_mesh("Armadillo_post_wmark.bin", m.verts, watermark_packed_inds);

    return 0;
}