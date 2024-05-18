<script>
  export let name = "";
  // export let value = document.getElementById("{name}-entry")[0].value;
  // {}
</script>

<div class="preamble">
  <span>
    {name} = 
  </span>
  <input id="{name}-entry" type="number" value="1">
</div>

<style>
  input {
    width:40px;
  }
  .preamble {
    display: flex;
    flex-flow: row;
    justify-content: space-between;
    width: 74px;
  }
</style>