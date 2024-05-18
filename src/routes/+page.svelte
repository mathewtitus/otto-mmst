<script>
  import Entry from "./Entry.svelte";
  import { onMount } from 'svelte';
  
  let p, q, r;
  // let setVars;

  onMount(() => {
    function setVars() {
      p = document.getElementById("p-entry").value;
      q = document.getElementById("q-entry").value;
      r = document.getElementById("r-entry").value;
    };

    setVars();
    document.addEventListener("input", setVars);
  })

</script>

<h1>Mean Minimal Spanning Tree Calculator</h1>

<row style="display: flex; flex-flow: row; justify-content: space-around;">
<column>
<Entry name="p"/>
<Entry name="q"/>
<Entry name="r"/>
</column>
<div>V = </div>
<column>
0 {p} {q} {r}<br>{p} 0 {r} {q}<br>{q} {r} 0 {p}<br>{r} {q} {p} 0
</column>
</row>

